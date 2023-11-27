version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf
import "../utilities/wf_metaspades_assembly.wdl" as metaspades_assembly_wf
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken_task
import "../../tasks/taxon_id/task_krona.wdl" as krona_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiameta_illumina_pe {
  meta {
    description: "Reference-based consensus calling or de novo assembly for metagenomic sequencing data"
  }
  input {
    File read1
    File read2
    String samplename
    File? reference
    File kraken2_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/k2_standard_08gb_20230605.tar.gz"
    Boolean output_additional_files = false
  }
  call kraken_task.kraken2_standalone as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = kraken2_db,
      kraken2_args = "",
      classified_out = "classified#.fastq",
      unclassified_out = "unclassified#.fastq"
  }
  call krona_task.krona as krona_raw {
    input:
      kraken2_report = kraken2_raw.kraken2_report,
      samplename = samplename
  }
  call read_qc_wf.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1_raw = read1,
        read2_raw = read2,
        workflow_series = "theiameta"
    }
  call kraken_task.kraken2_standalone as kraken2_clean {
    input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      kraken2_db = kraken2_db,
      kraken2_args = "",
      classified_out = "classified#.fastq",
      unclassified_out = "unclassified#.fastq"
  }
  call krona_task.krona as krona_clean {
    input:
      kraken2_report = kraken2_clean.kraken2_report,
      samplename = samplename
  }
  call metaspades_assembly_wf.metaspades_assembly_pe as metaspades {
    input:
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      samplename = samplename
    }
    # if reference is provided, perform mapping of assembled contigs to 
    # reference with minimap2, and extract those as final assembly
    if (defined(reference)){
      call minimap2_task.minimap2 as minimap2_assembly {
        input:
          query1 = metaspades.assembly_fasta,
          reference = select_first([reference]),
          samplename = samplename,
          mode = "asm20",
          output_sam = false
      }
      call parse_mapping_task.retrieve_aligned_contig_paf {
        input:
          paf = minimap2_assembly.minimap2_out,
          assembly = metaspades.assembly_fasta,
          samplename = samplename
      }
      call parse_mapping_task.calculate_coverage_paf {
        input:
          paf = minimap2_assembly.minimap2_out
      }
    }
    call quast_task.quast {
      input:
        assembly = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta]),
        samplename = samplename,
        min_contig_len = 1
      }
    if (output_additional_files){
      call minimap2_task.minimap2 as minimap2_reads {
        input:
          query1 = read_QC_trim.read1_clean,
          query2 = read_QC_trim.read2_clean, 
          reference = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta]),
          samplename = samplename,
          mode = "sr",
          output_sam = true
      }
      call parse_mapping_task.sam_to_sorted_bam {
        input:
          sam = minimap2_reads.minimap2_out,
          samplename = samplename
      }
      call parse_mapping_task.calculate_coverage {
        input:
          bam = sam_to_sorted_bam.bam,
          bai = sam_to_sorted_bam.bai
      }
      call parse_mapping_task.retrieve_pe_reads_bam as retrieve_unaligned_pe_reads_sam {
        input:
          bam = sam_to_sorted_bam.bam,
          samplename = samplename,
          prefix = "unassembled",
          sam_flag = 4
      }
      call parse_mapping_task.retrieve_pe_reads_bam as retrieve_aligned_pe_reads_sam {
        input:
          bam = sam_to_sorted_bam.bam,
          samplename = samplename,
          sam_flag = 2,
          prefix = "assembled"
      }
      call parse_mapping_task.assembled_reads_percent {
        input:
          bam = sam_to_sorted_bam.bam,
      } 
    }
    call versioning.version_capture{
      input:
  }
  output {
    # Version capture
    String theiameta_illumina_pe_version = version_capture.phb_version
    String theiameta_illumina_pe_analysis_date = version_capture.date
    # Kraken2 outputs
    String kraken2_version = kraken2_raw.kraken2_version
    String kraken2_docker = kraken2_raw.kraken2_docker
    File kraken2_report_raw = kraken2_raw.kraken2_report
    Float kraken2_percent_human_raw = kraken2_raw.kraken2_percent_human
    File kraken2_report_clean = kraken2_clean.kraken2_report
    Float kraken2_percent_human_clean = kraken2_clean.kraken2_percent_human
    # Krona outputs
    String krona_version = krona_raw.krona_version
    String krona_docker = krona_raw.krona_docker
    File krona_html_raw = krona_raw.krona_html
    File krona_html_clean = krona_clean.krona_html
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # Read QC - fastq_scan outputs
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    String fastq_scan_docker = read_QC_trim.fastq_scan_docker
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - bbduk outputs
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - Read stats
    Float? average_read_length = read_QC_trim.average_read_length
    # Assembly - metaspades 
    File assembly_fasta = select_first([retrieve_aligned_contig_paf.final_assembly, metaspades.assembly_fasta])
    String metaspades_version = metaspades.metaspades_version
    String metaspades_docker = metaspades.metaspades_docker
    # Assembly - minimap2
    String minimap2_version = metaspades.minimap2_version
    String minimap2_docker = metaspades.minimap2_docker
    # Assembly - samtools
    String samtools_version = metaspades.samtools_version
    String samtools_docker = metaspades.samtools_docker
    # Assembly - pilon
    String pilon_version = metaspades.pilon_version
    String pilon_docker = metaspades.pilon_docker
    # Assembly QC - quast
    Int assembly_length = quast.genome_length
    Int contig_number = quast.number_contigs
    Int largest_contig = quast.largest_contig
    String quast_version = quast.version
    String quast_docker = quast.quast_docker
    # Assembly QC - minimap2
    Float? percent_coverage = calculate_coverage_paf.percent_coverage
    # Assembly QC - bedtools
    Float? assembly_mean_coverage = calculate_coverage.mean_depth_coverage
    String? bedtools_version = calculate_coverage.bedtools_version
    String? bedtools_docker = calculate_coverage.bedtools_docker
    # Read retrieval
    File? read1_unmapped = retrieve_unaligned_pe_reads_sam.read1
    File? read2_unmapped = retrieve_unaligned_pe_reads_sam.read2
    File? read1_mapped = retrieve_aligned_pe_reads_sam.read1 
    File? read2_mapped = retrieve_aligned_pe_reads_sam.read2
    # Assembly stats
    Float? percentage_mapped_reads = assembled_reads_percent.percentage_mapped
    }
}