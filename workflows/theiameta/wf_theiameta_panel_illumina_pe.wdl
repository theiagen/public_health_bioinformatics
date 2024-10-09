version 1.0

import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/assembly/task_metaspades.wdl" as metaspades_task
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/read_filtering/task_pilon.wdl" as pilon_task
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken_task
import "../../tasks/taxon_id/task_krakentools.wdl" as krakentools_task
import "../../tasks/taxon_id/contamination/task_krona.wdl" as krona_task
import "../../tasks/utilities/data_handling/task_gather_scatter.wdl" as gather_scatter_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_workflow
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_trim_pe

workflow theiameta_panel_illumina_pe {
  input {
    String samplename
    File read1
    File read2
    Array[Int]? taxon_ids # suggest using a workspace element if user wants to modify?

    Int minimum_read_number = 100
    File kraken2_db = "gs://theiagen-large-public-files-rp/terra/databases/kraken2/k2_viral_20240112.tar.gz"
  }
  # kraken does not run as part of the theiameta track in read_QC_trim -- we may want to change that
  call kraken_task.kraken2_standalone as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = kraken2_db
  }
  call krona_task.krona as krona_raw {
    input:
      kraken2_report = kraken2_raw.kraken2_report,
      samplename = samplename
  }
  call read_qc_trim_pe.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2,
        workflow_series = "theiameta"
  }
  call kraken_task.kraken2_standalone as kraken2_clean {
    input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      kraken2_db = kraken2_db
  }
  call krona_task.krona as krona_clean {
    input:
      kraken2_report = kraken2_clean.kraken2_report,
      samplename = samplename
  }
  scatter (taxon_id in taxon_ids) {
    call krakentools_task.extract_kraken_reads as krakentools {
      input:
        # we should consider changing the classified_report name so 
        #  it won't be confused with the actual kraken2 report
        kraken2_output = kraken2_clean.kraken2_classified_report,
        kraken2_report = kraken2_clean.kraken2_report,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        taxon_id = taxon_id
    }
    if (krakentools.success) {
      call fastq_scan.fastq_scan_pe as fastq_scan_binned {
        input:
          read1 = krakentools.extracted_read1,
          read2 = krakentools.extracted_read2
      }
      #### ADJUST IN THE FUTURE; SETTING TO 100 FOR TESTING ####
      if (fastq_scan_binned.read1_seq > ~{minimum_read_number}) {
        String did_attempt_assembly = "Assembly attempted"
        call metaspades_task.metaspades_pe {
          input:
            read1_cleaned = krakentools.extracted_read1,
            read2_cleaned = krakentools.extracted_read2,
            samplename = "~{samplename}_~{taxon_id}"
        }
        call minimap2_task.minimap2 as minimap2_assembly_correction {
          input:
            query1 = krakentools.extracted_read1,
            query2 = krakentools.extracted_read2,
            reference = metaspades_pe.assembly_fasta,
            samplename = "~{samplename}_~{taxon_id}",
            mode = "sr",
            output_sam = true
        }
        call parse_mapping_task.sam_to_sorted_bam as sort_bam_assembly_correction {
          input:
            sam = minimap2_assembly_correction.minimap2_out,
            samplename = "~{samplename}_~{taxon_id}"
        }
        call pilon_task.pilon {
          input:
            assembly = metaspades_pe.assembly_fasta,
            bam = sort_bam_assembly_correction.bam,
            bai = sort_bam_assembly_correction.bai,
            samplename = "~{samplename}_~{taxon_id}"
        }    
        call quast_task.quast {
          input:
            assembly = pilon.assembly_fasta,
            samplename = "~{samplename}_~{taxon_id}",
            min_contig_length = 1
        }
        call morgana_magic_workflow.morgana_magic {
          input:
            samplename = "~{samplename}_~{taxon_id}",
            assembly_fasta = pilon.assembly_fasta,
            read1 = krakentools.extracted_read1,
            read2 = krakentools.extracted_read2,
            taxon_id = "~{taxon_id}"
        }
      }
    }
  }  
  call gather_scatter_task.gather_scatter {
    input:
      taxon_ids = select_first([taxon_ids]),
      organism = morgana_magic.organism,
      extracted_read1 = krakentools.extracted_read1,
      extracted_read2 = krakentools.extracted_read2,
      fastq_scan_num_reads_binned1 = fastq_scan_binned.read1_seq,
      fastq_scan_num_reads_binned2 = fastq_scan_binned.read2_seq,
      fastq_scan_num_reads_binned_pairs = fastq_scan_binned.read_pairs,
      pilon_assembly_fasta = pilon.assembly_fasta,
      quast_genome_length = quast.genome_length,
      quast_number_contigs = quast.number_contigs,
      quast_n50 = quast.n50_value,
      quast_gc_percent = quast.gc_percent,
      pango_lineage = morgana_magic.pango_lineage,
      pango_lineage_expanded = morgana_magic.pango_lineage_expanded,
      pangolin_conflicts = morgana_magic.pangolin_conflicts,
      pangolin_notes = morgana_magic.pangolin_notes,
      pangolin_assignment_version = morgana_magic.pangolin_assignment_version,
      pangolin_versions = morgana_magic.pangolin_versions,
      pangolin_docker = morgana_magic.pangolin_docker,
      nextclade_version = morgana_magic.nextclade_version,
      nextclade_docker = morgana_magic.nextclade_docker,
      nextclade_ds_tag = morgana_magic.nextclade_ds_tag,
      nextclade_aa_subs = morgana_magic.nextclade_aa_subs,
      nextclade_aa_dels = morgana_magic.nextclade_aa_dels,
      nextclade_clade = morgana_magic.nextclade_clade,
      nextclade_lineage = morgana_magic.nextclade_lineage,
      nextclade_qc = morgana_magic.nextclade_qc
  } 
  output {
    Array[String?] identified_organisms = select_first([morgana_magic.organism])
  }
}