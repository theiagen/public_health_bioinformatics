version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/quality_control/advanced_metrics/task_checkv.wdl" as checkv_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/taxon_id/task_skani.wdl" as skani_task
import "../../tasks/utilities/data_import/task_ncbi_datasets.wdl" as ncbi_datasets_task
import "../../tasks/taxon_id/task_identify_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as ivar_consensus
import "../../tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf
import "../../tasks/taxon_id/task_krakentools.wdl" as krakentools_task
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2_task
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/utilities/file_handling/task_cat_lanes.wdl" as cat_lanes
import "../../tasks/utilities/data_export/task_export_taxon_table.wdl" as export_taxon_table_task

workflow theiaviral_panel {
  input {
    File read1
    File read2
    String samplename
    Array[String] taxon_ids
    File output_taxon_table

    String terra_project
    String terra_workspace
    String kraken_db

    Int min_map_quality = 20
    Int min_depth = 10
    Float min_allele_freq = 0.6


    Boolean extract_unclassified = false
    Int minimum_read_number = 1000
    Boolean skip_metaviralspades = false

  }  
  call versioning.version_capture {
    input:
  }
  # read QC, classification, extraction, and trimming
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      kraken_db = kraken_db,
      workflow_series = "theiacov" # this will return kraken2 reports that are necessary
  }
  call kraken2_task.kraken2_standalone as kraken2 {
    input:
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      samplename = samplename,
      kraken2_db = kraken_db
  }
  # get kraken outputs ready will have to run it outside of read qc trim to fix
  scatter (taxon_id in taxon_ids) {
    call krakentools_task.extract_kraken_reads as krakentools {
      input:
        kraken2_output = kraken2.kraken2_classified_report,
        kraken2_report = kraken2.kraken2_report,
        read1 = kraken2.kraken2_classified_read1,
        read2 = kraken2.kraken2_classified_read2,
        taxon_id = taxon_id
    }
    if (krakentools.success) {
      if (extract_unclassified) {
        call cat_lanes.cat_lanes {
          input:
            samplename = samplename + "_" + taxon_id,
            read1_lane1 = kraken2.kraken2_unclassified_read1,
            read1_lane2 = select_first([krakentools.extracted_read1]),
            read2_lane1 = kraken2.kraken2_unclassified_read2,
            read2_lane2 = select_first([krakentools.extracted_read2])
        }
      }
      call fastq_scan.fastq_scan_pe as fastq_scan_binned {
        input:
          read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
          read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2])
      }
      if (fastq_scan_binned.read1_seq > minimum_read_number) {
        call ncbi_datasets_task.ncbi_datasets_viral_taxon_summary as ncbi_taxon_summary {
          input:
            taxon = taxon_id
        }
        call identify_taxon_id_task.identify_taxon_id as ncbi_identify {
          input:
            taxon = taxon_id
        }
        if (! skip_metaviralspades) {
          call spades_task.spades {
            input:
              read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
              samplename =  samplename + "_" + taxon_id,
              spades_type = "metaviral"
          }
        }
        # fallback to megahit if metaviralspades fails to identify a complete virus
        if (select_first([spades.spades_status, "FAIL"]) == "FAIL") {
          call megahit_task.megahit {
            input:
              read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
              samplename =  samplename + "_" + taxon_id,
              hard_fail = false
          }
        }
        # check if assembly exists
        if (defined(spades.assembly_fasta) || defined(megahit.assembly_fasta)) {
          # quality control metrics for de novo assembly (ie. completeness, viral gene count, contamination)
          call checkv_task.checkv as checkv_denovo {
            input:
              assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
              samplename =  samplename + "_" + taxon_id
          }
          # quality control metrics for de novo assembly (ie. contigs, n50, GC content, genome length)
          call quast_task.quast as quast_denovo {
            input:
              assembly = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
              samplename =  samplename + "_" + taxon_id
          }
          # ANI-based reference genome selection
          # do we want a skani failure to be a hard failure? or be like previous kraken failures?
          call skani_task.skani as skani {
            input:
              assembly_fasta = select_first([spades.assembly_fasta, megahit.assembly_fasta]),
              samplename =  samplename + "_" + taxon_id
          }
          # download the best reference determined from skani
          call ncbi_datasets_task.ncbi_datasets_download_genome_accession as ncbi_datasets {
            input:
              ncbi_accession = skani.skani_top_accession,
              use_ncbi_virus = true
          }
          # align reads to reference
          call bwa_task.bwa {
            input:
              samplename =  samplename + "_" + taxon_id,
              read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              read2 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              reference_genome = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta])
          }
          # consensus calling via ivar
          call ivar_consensus.consensus {
            input:
              bamfile = bwa.sorted_bam,
              samplename =  samplename + "_" + taxon_id,
              reference_genome = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta]),
              min_qual = min_map_quality,
              consensus_min_depth = select_first([min_depth, 10]),
              consensus_min_freq = min_allele_freq,
              all_positions = true
          }
          # variant calling via ivar
          call variant_call_task.variant_call as ivar_variants {
            input:
              mpileup = consensus.sample_mpileup,
              samplename =  samplename + "_" + taxon_id,
              reference_genome = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta]),
              min_qual = min_map_quality,
              organism = "",
              variant_min_freq = min_allele_freq,
              variant_min_depth = select_first([min_depth, 10])
          }
          # quality control metrics for reads mapping to reference (ie. coverage, depth, base/map quality)
          call assembly_metrics_task.stats_n_coverage as read_mapping_stats {
            input:
              bamfile = bwa.sorted_bam,
              samplename =  samplename + "_" + taxon_id
          }
          # quality control metrics for consensus (ie. number of bases, degenerate bases, genome length)
          call consensus_qc_task.consensus_qc as consensus_qc {
            input:
              assembly_fasta = consensus.consensus_seq,
              reference_genome = select_first([ncbi_datasets.ncbi_datasets_assembly_fasta]),
              genome_length = select_first([ncbi_taxon_summary.avg_genome_length])
          }
          # quality control metrics for consensus (ie. completeness, viral gene count, contamination)
          call checkv_task.checkv as checkv_consensus {
            input:
              assembly = consensus.consensus_seq,
              samplename =  samplename + "_" + taxon_id
          }
          # morgana magic
          call morgana_magic_wf.morgana_magic {
            input:
              samplename = samplename + "_" + taxon_id,
              assembly_fasta = consensus.consensus_seq,
              read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              read2 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              taxon_name = ncbi_identify.taxon_name,
              seq_method = "panel"
          }
        }
        # call export_taxon_table
        call export_taxon_table_task.export_taxon_table {
          input:
            samplename = samplename + "_" + taxon_id,
            taxon_table = output_taxon_table,
            gambit_predicted_taxon = ncbi_identify.taxon_name,
            terra_project = terra_project,
            terra_workspace = terra_workspace,
            columns_to_export = {
              "samplename": samplename + "_" + taxon_id,
              "taxon_name": ncbi_identify.taxon_name,
              "assembly_fasta": consensus.consensus_seq,
              "read1": select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              "read2": select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              "fastq_scan_num_reads_binned1": fastq_scan_binned.read1_seq,
              "fastq_scan_num_reads_binned2": fastq_scan_binned.read2_seq,
              "fastq_scan_num_reads_binned_pairs": fastq_scan_binned.read_pairs,
              "fastq_scan_docker": fastq_scan_binned.fastq_scan_docker,
              "fastq_scan_version": fastq_scan_binned.version,
              "fastq_scan_binned1_json": fastq_scan_binned.read1_fastq_scan_json,
              "fastq_scan_binned2_json": fastq_scan_binned.read2_fastq_scan_json,
              "ncbi_taxon_summary_tsv": ncbi_taxon_summary.taxon_summary_tsv,
              "ncbi_taxon_summary_avg_genome_length": ncbi_taxon_summary.avg_genome_length,
              "assembly_denovo_fasta": select_first([spades.assembly_fasta, megahit.assembly_fasta, ""]),
              "metaviralspades_status": spades.spades_status,
              "metaviralspades_version": spades.spades_version,
              "metaviralspades_docker": spades.spades_docker,
              "megahit_version": megahit.megahit_version,
              "megahit_docker": megahit.megahit_docker,
              "checkv_denovo_summary": checkv_denovo.checkv_summary,
              "checkv_denovo_contamination": checkv_denovo.checkv_contamination,
              "checkv_denovo_weighted_completeness": checkv_denovo.weighted_completeness,
              "checkv_denovo_weighted_contamination": checkv_denovo.weighted_contamination,
              "checkv_denovo_total_genes": checkv_denovo.total_genes,
              "checkv_denovo_version": checkv_denovo.checkv_version,
              "quast_denovo_report": quast_denovo.quast_report,
              "quast_denovo_genome_length": quast_denovo.genome_length,
              "quast_denovo_number_contigs": quast_denovo.number_contigs,
              "quast_denovo_n50_value": quast_denovo.n50_value,
              "quast_denovo_largest_contig": quast_denovo.largest_contig,
              "quast_denovo_gc_percent": quast_denovo.gc_percent,
              "quast_denovo_uncalled_bases": quast_denovo.uncalled_bases,
              "quast_denovo_version": quast_denovo.version,
              "quast_denovo_docker": quast_denovo.quast_docker,
              "skani_report": skani.skani_report,
              "skani_top_accession": skani.skani_top_accession,
              "skani_top_score": skani.skani_top_score,
              "skani_top_ani": skani.skani_top_ani,
              "skani_top_ref_coverage": skani.skani_top_ref_coverage,
              "skani_database": skani.skani_database,
              "skani_version": skani.skani_version,
              "skani_docker": skani.skani_docker,
              "skani_top_ani_fasta": ncbi_datasets.ncbi_datasets_assembly_fasta,
              "bwa_version": bwa.bwa_version,
              "bwa_samtools_version": bwa.sam_version,
              "bwa_read1_aligned": bwa.read1_aligned,
              "bwa_read2_aligned": bwa.read2_aligned,
              "bwa_sorted_bam": bwa.sorted_bam,
              "bwa_sorted_bai": bwa.sorted_bai,
              "bwa_read1_unaligned": bwa.read1_unaligned,
              "bwa_read2_unaligned": bwa.read2_unaligned,
              "bwa_sorted_bam_unaligned": bwa.sorted_bam_unaligned,
              "bwa_sorted_bai_unaligned": bwa.sorted_bam_unaligned_bai,
              "ivar_tsv": ivar_variants.sample_variants_tsv,
              "ivar_vcf": ivar_variants.sample_variants_vcf,
              "ivar_variant_version": ivar_variants.ivar_version,
              "read_mapping_report": read_mapping_stats.metrics_txt,
              "read_mapping_statistics": read_mapping_stats.stats,
              "read_mapping_cov_hist": read_mapping_stats.cov_hist,
              "read_mapping_cov_stats": read_mapping_stats.cov_stats,
              "read_mapping_flagstat": read_mapping_stats.flagstat,
              "read_mapping_coverage": read_mapping_stats.coverage,
              "read_mapping_depth": read_mapping_stats.depth,
              "read_mapping_meanbaseq": read_mapping_stats.meanbaseq,
              "read_mapping_meanmapq": read_mapping_stats.meanmapq,
              "read_mapping_percentage_mapped_reads": read_mapping_stats.percentage_mapped_reads,
              "read_mapping_samtools_version": read_mapping_stats.samtools_version,
              "consensus_qc_number_N": consensus_qc.number_N,
              "consensus_qc_assembly_length_unambiguous": consensus_qc.number_ATCG,
              "consensus_qc_number_degenerate": consensus_qc.number_Degenerate,
              "consensus_qc_number_total": consensus_qc.number_Total,
              "consensus_qc_percent_reference_coverage": consensus_qc.percent_reference_coverage,
              "checkv_consensus_summary": checkv_consensus.checkv_summary,
              "checkv_consensus_contamination": checkv_consensus.checkv_contamination,
              "checkv_consensus_weighted_completeness": checkv_consensus.weighted_completeness,
              "checkv_consensus_weighted_contamination": checkv_consensus.weighted_contamination,
              "checkv_consensus_total_genes": checkv_consensus.total_genes,
              "checkv_consensus_version": checkv_consensus.checkv_version,
              "organism": morgana_magic.organism,
              "pango_lineage": morgana_magic.pango_lineage,
              "pango_lineage_expanded": morgana_magic.pango_lineage_expanded,
              "pangolin_conflicts": morgana_magic.pangolin_conflicts,
              "pangolin_notes": morgana_magic.pangolin_notes,
              "pangolin_assignment_version": morgana_magic.pangolin_assignment_version,
              "pango_lineage_report": morgana_magic.pango_lineage_report,
              "pangolin_docker": morgana_magic.pangolin_docker,
              "pangolin_versions": morgana_magic.pangolin_versions,
              "nextclade_version": morgana_magic.nextclade_version,
              "nextclade_docker": morgana_magic.nextclade_docker,
              "nextclade_json": morgana_magic.nextclade_json,
              "auspice_json": morgana_magic.auspice_json,
              "nextclade_tsv": morgana_magic.nextclade_tsv,
              "nextclade_ds_tag": morgana_magic.nextclade_ds_tag,
              "nextclade_aa_subs": morgana_magic.nextclade_aa_subs,
              "nextclade_aa_dels": morgana_magic.nextclade_aa_dels,
              "nextclade_clade": morgana_magic.nextclade_clade,
              "nextclade_lineage": morgana_magic.nextclade_lineage,
              "nextclade_qc": morgana_magic.nextclade_qc,
              "nextclade_json_flu_ha": morgana_magic.nextclade_json_flu_ha,
              "auspice_json_flu_ha": morgana_magic.auspice_json_flu_ha,
              "nextclade_tsv_flu_ha": morgana_magic.nextclade_tsv_flu_ha,
              "nextclade_ds_tag_flu_ha": morgana_magic.nextclade_ds_tag_flu_ha,
              "nextclade_aa_subs_flu_ha": morgana_magic.nextclade_aa_subs_flu_ha,
              "nextclade_aa_dels_flu_ha": morgana_magic.nextclade_aa_dels_flu_ha,
              "nextclade_clade_flu_ha": morgana_magic.nextclade_clade_flu_ha,
              "nextclade_qc_flu_ha": morgana_magic.nextclade_qc_flu_ha,
              "nextclade_json_flu_na": morgana_magic.nextclade_json_flu_na,
              "auspice_json_flu_na": morgana_magic.auspice_json_flu_na,
              "nextclade_tsv_flu_na": morgana_magic.nextclade_tsv_flu_na,
              "nextclade_ds_tag_flu_na": morgana_magic.nextclade_ds_tag_flu_na,
              "nextclade_aa_subs_flu_na": morgana_magic.nextclade_aa_subs_flu_na,
              "nextclade_aa_dels_flu_na": morgana_magic.nextclade_aa_dels_flu_na,
              "nextclade_clade_flu_na": morgana_magic.nextclade_clade_flu_na,
              "nextclade_qc_flu_na": morgana_magic.nextclade_qc_flu_na,
              "irma_version": morgana_magic.irma_version,
              "irma_docker": morgana_magic.irma_docker,
              "irma_type": morgana_magic.irma_type,
              "irma_subtype": morgana_magic.irma_subtype,
              "irma_subtype_notes": morgana_magic.irma_subtype_notes,
              "genoflu_version": morgana_magic.genoflu_version,
              "genoflu_genotype": morgana_magic.genoflu_genotype,
              "genoflu_all_segments": morgana_magic.genoflu_all_segments,
              "genoflu_output_tsv": morgana_magic.genoflu_output_tsv,
              "abricate_flu_type": morgana_magic.abricate_flu_type,
              "abricate_flu_subtype":  morgana_magic.abricate_flu_subtype,
              "abricate_flu_results": morgana_magic.abricate_flu_results,
              "abricate_flu_database":  morgana_magic.abricate_flu_database,
              "abricate_flu_version": morgana_magic.abricate_flu_version
            }
        }
      }
    }  
  }
 
  output {
    # Number of assembled viruses
    Int assembled_viruses = length(select_all(export_taxon_table.status))
    # Workflow Versioning
    String theiaviral_panel_version = version_capture.phb_version
    String theiaviral_pannel_analysis_date = version_capture.date
    # Standalone Kraken2 outputs
    String kraken2_version = kraken2.kraken2_version
    String kraken2_database = kraken2.kraken2_database
    String kraken2_docker = kraken2.kraken2_docker
    File kraken2_report = kraken2.kraken2_report
    File kraken2_classified_report = kraken2.kraken2_classified_report
  }
}