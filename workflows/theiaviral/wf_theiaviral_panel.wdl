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
import "wf_theiaviral_illumina_pe.wdl" as theiaviral_illumina_pe

workflow theiaviral_panel {
  input {
    File read1
    File read2
    String samplename
    Array[String] taxon_ids
    File output_taxon_table

    String terra_project
    String terra_workspace
    String kraken_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/k2_viral_20240112.tar.gz"
    File? skani_db
    File? checkv_db

    Int min_map_quality = 20
    Int min_depth = 10
    Float min_allele_freq = 0.6
    String? read_extraction_rank

    Boolean extract_unclassified = false
    Boolean skip_theiaviral_screen = false
    Int minimum_read_number = 1000
    Boolean call_metaviralspades = false
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
        read2 = select_first([kraken2.kraken2_classified_read2]),
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
        # get the taxon information from ncbi
        call identify_taxon_id_task.identify_taxon_id as ncbi_identify {
          input:
            taxon = taxon_id,
            rank = read_extraction_rank,
            use_ncbi_virus = true
        }
        call theiaviral_illumina_pe.theiaviral_illumina_pe as theiaviral {
          input:
            read1 = select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
            read2 = select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
            samplename = samplename + "_" + taxon_id,
            taxon = taxon_id,
            call_metaviralspades = call_metaviralspades,
            kraken_db = kraken_db,
            skip_qc = true,
            skip_screen = skip_theiaviral_screen,
            skani_db = skani_db,
            checkv_db = checkv_db,
            genome_length = ncbi_identify.avg_genome_length,
            min_map_quality = min_map_quality,
            min_depth = min_depth,
            min_allele_freq = min_allele_freq,
            read_extraction_rank = read_extraction_rank
        }
        # call export_taxon_table
        call export_taxon_table_task.export_taxon_table_vsp {
          input:
            samplename = samplename + "_" + sub(ncbi_identify.taxon_name, " ", "-"),
            taxon_table = output_taxon_table,
            gambit_predicted_taxon = ncbi_identify.taxon_name,
            terra_project = terra_project,
            terra_workspace = terra_workspace,
            columns_to_export = {
              "samplename": samplename + "_" + taxon_id,
              "taxon_name": ncbi_identify.taxon_name,
              "assembly_fasta": theiaviral.assembly_consensus_fasta,
              "read1": select_first([cat_lanes.read1_concatenated, krakentools.extracted_read1]),
              "read2": select_first([cat_lanes.read2_concatenated, krakentools.extracted_read2]),
              "fastq_scan_num_reads_binned1": fastq_scan_binned.read1_seq,
              "fastq_scan_num_reads_binned2": fastq_scan_binned.read2_seq,
              "fastq_scan_num_reads_binned_pairs": fastq_scan_binned.read_pairs,
              "fastq_scan_docker": fastq_scan_binned.fastq_scan_docker,
              "fastq_scan_version": fastq_scan_binned.version,
              "fastq_scan_binned1_json": fastq_scan_binned.read1_fastq_scan_json,
              "fastq_scan_binned2_json": fastq_scan_binned.read2_fastq_scan_json,
              "ncbi_taxon_summary_tsv": ncbi_identify.taxon_summary_tsv,
              "ncbi_taxon_summary_avg_genome_length": ncbi_identify.avg_genome_length,
              "assembly_denovo_fasta": theiaviral.assembly_denovo_fasta,
              "metaviralspades_status": theiaviral.metaviralspades_status,
              "metaviralspades_version": theiaviral.metaviralspades_version,
              "metaviralspades_docker": theiaviral.metaviralspades_docker,
              "megahit_version": theiaviral.megahit_version,
              "megahit_docker": theiaviral.megahit_docker,
              "checkv_denovo_summary": theiaviral.checkv_denovo_summary,
              "checkv_denovo_contamination": theiaviral.checkv_denovo_contamination,
              "checkv_denovo_weighted_completeness": theiaviral.checkv_denovo_weighted_completeness,
              "checkv_denovo_weighted_contamination": theiaviral.checkv_denovo_weighted_contamination,
              "checkv_denovo_total_genes": theiaviral.checkv_denovo_total_genes,
              "checkv_denovo_version": theiaviral.checkv_denovo_version,
              "quast_denovo_report": theiaviral.quast_denovo_report,
              "quast_denovo_genome_length": theiaviral.quast_denovo_genome_length,
              "quast_denovo_number_contigs": theiaviral.quast_denovo_number_contigs,
              "quast_denovo_n50_value": theiaviral.quast_denovo_n50_value,
              "quast_denovo_largest_contig": theiaviral.quast_denovo_largest_contig,
              "quast_denovo_gc_percent": theiaviral.quast_denovo_gc_percent,
              "quast_denovo_uncalled_bases": theiaviral.quast_denovo_uncalled_bases,
              "quast_denovo_version": theiaviral.quast_denovo_version,
              "quast_denovo_docker": theiaviral.quast_denovo_docker,
              "skani_report": theiaviral.skani_report,
              "skani_top_accession": theiaviral.skani_top_accession,
              "skani_top_score": theiaviral.skani_top_score,
              "skani_top_ani": theiaviral.skani_top_ani,
              "skani_database": theiaviral.skani_database,
              "skani_version": theiaviral.skani_version,
              "skani_docker": theiaviral.skani_docker,
              "skani_top_ani_fasta": theiaviral.skani_top_ani_fasta,
              "bwa_version": theiaviral.bwa_version,
              "bwa_samtools_version": theiaviral.bwa_samtools_version,
              "bwa_read1_aligned": theiaviral.bwa_read1_aligned,
              "bwa_read2_aligned": theiaviral.bwa_read2_aligned,
              "bwa_sorted_bam": theiaviral.bwa_sorted_bam,
              "bwa_sorted_bai": theiaviral.bwa_sorted_bai,
              "bwa_read1_unaligned": theiaviral.bwa_read1_unaligned,
              "bwa_read2_unaligned": theiaviral.bwa_read2_unaligned,
              "bwa_sorted_bam_unaligned": theiaviral.bwa_sorted_bam_unaligned,
              "bwa_sorted_bai_unaligned": theiaviral.bwa_sorted_bam_unaligned_bai,
              "ivar_tsv": theiaviral.ivar_tsv,
              "ivar_vcf": theiaviral.ivar_vcf,
              "ivar_variant_version": theiaviral.ivar_variant_version,
              "read_mapping_report": theiaviral.read_mapping_report,
              "read_mapping_statistics": theiaviral.read_mapping_statistics,
              "read_mapping_cov_hist": theiaviral.read_mapping_cov_hist,
              "read_mapping_cov_stats": theiaviral.read_mapping_cov_stats,
              "read_mapping_flagstat": theiaviral.read_mapping_flagstat,
              "read_mapping_coverage": theiaviral.read_mapping_coverage,
              "read_mapping_depth": theiaviral.read_mapping_depth,
              "read_mapping_meanbaseq": theiaviral.read_mapping_meanbaseq,
              "read_mapping_meanmapq": theiaviral.read_mapping_meanmapq,
              "read_mapping_percentage_mapped_reads": theiaviral.read_mapping_percentage_mapped_reads,
              "read_mapping_samtools_version": theiaviral.read_mapping_samtools_version,
              "consensus_qc_number_N": theiaviral.consensus_qc_number_N,
              "consensus_qc_assembly_length_unambiguous": theiaviral.consensus_qc_assembly_length_unambiguous,
              "consensus_qc_number_degenerate": theiaviral.consensus_qc_number_Degenerate,
              "consensus_qc_number_total": theiaviral.consensus_qc_number_Total,
              "consensus_qc_percent_reference_coverage": theiaviral.consensus_qc_percent_reference_coverage,
              "checkv_consensus_summary": theiaviral.checkv_consensus_summary,
              "checkv_consensus_contamination": theiaviral.checkv_consensus_contamination,
              "checkv_consensus_weighted_completeness": theiaviral.checkv_consensus_weighted_completeness,
              "checkv_consensus_weighted_contamination": theiaviral.checkv_consensus_weighted_contamination,
              "checkv_consensus_total_genes": theiaviral.checkv_consensus_total_genes,
              "checkv_consensus_version": theiaviral.checkv_consensus_version,
              "pango_lineage": theiaviral.pango_lineage,
              "pango_lineage_expanded": theiaviral.pango_lineage_expanded,
              "pangolin_conflicts": theiaviral.pangolin_conflicts,
              "pangolin_notes": theiaviral.pangolin_notes,
              "pangolin_assignment_version": theiaviral.pangolin_assignment_version,
              "pango_lineage_report": theiaviral.pango_lineage_report,
              "pangolin_docker": theiaviral.pangolin_docker,
              "pangolin_versions": theiaviral.pangolin_versions,
              "nextclade_version": theiaviral.nextclade_version,
              "nextclade_docker": theiaviral.nextclade_docker,
              "nextclade_json_mpxv": theiaviral.nextclade_json_mpxv,
              "auspice_json_mpxv": theiaviral.auspice_json_mpxv,
              "nextclade_tsv_mpxv": theiaviral.nextclade_tsv_mpxv,
              "nextclade_ds_tag": theiaviral.nextclade_ds_tag,
              "nextclade_aa_subs_mpxv": theiaviral.nextclade_aa_subs_mpxv,
              "nextclade_aa_dels_mpxv": theiaviral.nextclade_aa_dels_mpxv,
              "nextclade_clade_mpxv": theiaviral.nextclade_clade_mpxv,
              "nextclade_lineage_mpxv": theiaviral.nextclade_lineage_mpxv,
              "nextclade_qc_mpxv": theiaviral.nextclade_qc_mpxv,
              "nextclade_json_flu_ha": theiaviral.nextclade_json_flu_ha,
              "auspice_json_flu_ha": theiaviral.auspice_json_flu_ha,
              "nextclade_tsv_flu_ha": theiaviral.nextclade_tsv_flu_ha,
              "nextclade_ds_tag_flu_ha": theiaviral.nextclade_ds_tag_flu_ha,
              "nextclade_aa_subs_flu_ha": theiaviral.nextclade_aa_subs_flu_ha,
              "nextclade_aa_dels_flu_ha": theiaviral.nextclade_aa_dels_flu_ha,
              "nextclade_clade_flu_ha": theiaviral.nextclade_clade_flu_ha,
              "nextclade_qc_flu_ha": theiaviral.nextclade_qc_flu_ha,
              "nextclade_json_flu_na": theiaviral.nextclade_json_flu_na,
              "auspice_json_flu_na": theiaviral.auspice_json_flu_na,
              "nextclade_tsv_flu_na": theiaviral.nextclade_tsv_flu_na,
              "nextclade_ds_tag_flu_na": theiaviral.nextclade_ds_tag_flu_na,
              "nextclade_aa_subs_flu_na": theiaviral.nextclade_aa_subs_flu_na,
              "nextclade_aa_dels_flu_na": theiaviral.nextclade_aa_dels_flu_na,
              "nextclade_clade_flu_na": theiaviral.nextclade_clade_flu_na,
              "nextclade_qc_flu_na": theiaviral.nextclade_qc_flu_na,
              "irma_version": theiaviral.irma_version,
              "irma_docker": theiaviral.irma_docker,
              "irma_type": theiaviral.irma_type,
              "irma_subtype": theiaviral.irma_subtype,
              "irma_subtype_notes": theiaviral.irma_subtype_notes,
              "genoflu_version": theiaviral.genoflu_version,
              "genoflu_genotype": theiaviral.genoflu_genotype,
              "genoflu_all_segments": theiaviral.genoflu_all_segments,
              "genoflu_output_tsv": theiaviral.genoflu_output_tsv,
              "abricate_flu_type": theiaviral.abricate_flu_type,
              "abricate_flu_subtype":  theiaviral.abricate_flu_subtype,
              "abricate_flu_results": theiaviral.abricate_flu_results,
              "abricate_flu_database":  theiaviral.abricate_flu_database,
              "abricate_flu_version": theiaviral.abricate_flu_version
            }
        }
      }
    }  
  }
 
  output {
    # Number of assembled viruses
    Int assembled_viruses = length(select_all(export_taxon_table_vsp.status))
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