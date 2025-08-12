version 1.0

import "../../tasks/assembly/task_artic_consensus.wdl" as artic_consensus
import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

workflow theiacov_clearlabs {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }
  input {
    String samplename
    File read1
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    # assembly parameters
    Int normalise = 20000
    String medaka_docker = "us-docker.pkg.dev/general-theiagen/staphb/artic-ncov2019:1.3.0-medaka-1.4.3"
    # reference values
    File? reference_genome
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # kraken parameters
    String? target_organism
    # qc check parameters
    File? qc_check_table
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      reference_genome = reference_genome,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,
      kraken_target_organism_input = target_organism
  }
  call fastq_scan.fastq_scan_se as fastq_scan_raw_reads {
    input:
      read1 = read1
  }
  call ncbi_scrub.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = read1
  }
  call fastq_scan.fastq_scan_se as fastq_scan_clean_reads {
    input:
      read1 = ncbi_scrub_se.read1_dehosted
  }
  call kraken2.kraken2_theiacov as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1,
      target_organism = organism_parameters.kraken_target_organism
  }  
  call kraken2.kraken2_theiacov as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted,
      target_organism = organism_parameters.kraken_target_organism
  }
  call artic_consensus.consensus {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted,
      primer_bed = primer_bed,
      normalise = normalise,
      docker = medaka_docker,
      reference_genome = organism_parameters.reference,
      organism = organism
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = consensus.sorted_bam
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = consensus.consensus_seq,
      reference_genome = organism_parameters.reference
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = consensus.consensus_seq,
        docker = organism_parameters.pangolin_docker
    }
    call gene_coverage_task.gene_coverage {
      input: 
        samplename = samplename,
        bamfile = consensus.trim_sorted_bam,
        bedfile = organism_parameters.gene_locations_bed,
        organism = organism_parameters.standardized_organism,
        min_depth = 20
    }
  }
  if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2") {
    # tasks specific to either MPXV or sars-cov-2
    call nextclade_task.nextclade_v3 {
      input:
      genome_fasta = consensus.consensus_seq,
      dataset_name = organism_parameters.nextclade_dataset_name,
      dataset_tag = organism_parameters.nextclade_dataset_tag
    }
    call nextclade_task.nextclade_output_parser {
      input:
      nextclade_tsv = nextclade_v3.nextclade_tsv,
      organism = organism
    }
    call vadr_task.vadr {
      input:
        genome_fasta = consensus.consensus_seq,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        max_length = organism_parameters.vadr_maxlength,
        vadr_opts = organism_parameters.vadr_opts,
        vadr_model_file = organism_parameters.vadr_model_file,
        skip_length = organism_parameters.vadr_skiplength,
        memory = organism_parameters.vadr_memory
    }
  }
  if (defined(qc_check_table)) {
    call qc_check.qc_check_phb as qc_check_task {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism_parameters.standardized_organism,
        num_reads_raw1 = fastq_scan_raw_reads.read1_seq,
        num_reads_clean1 = fastq_scan_clean_reads.read1_seq,
        kraken_human = kraken2_raw.percent_human,
        # kraken_sc2 = kraken2_raw.percent_sc2,
        # kraken_target_organism = kraken2_raw.percent_target_organism,
        kraken_human_dehosted = kraken2_dehosted.percent_human,
        # kraken_sc2_dehosted = kraken2_dehosted.percent_sc2,
        # kraken_target_organism_dehosted = kraken2_dehosted.percent_target_organism,
        meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq,
        assembly_mean_coverage = stats_n_coverage_primtrim.depth,
        number_N = consensus_qc.number_N,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        number_Degenerate =  consensus_qc.number_Degenerate,
        percent_reference_coverage =  consensus_qc.percent_reference_coverage,
        # sc2_s_gene_mean_coverage = gene_coverage.sc2_s_gene_depth,
        # sc2_s_gene_percent_coverage = gene_coverage.sc2_s_gene_percent_coverage,
        vadr_num_alerts = vadr.num_alerts
    }
  }  
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String theiacov_clearlabs_version = version_capture.phb_version
    String theiacov_clearlabs_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC - dehosting outputs
    File read1_dehosted = ncbi_scrub_se.read1_dehosted
    # Read QC - fastq_scan outputs
    Int fastq_scan_num_reads_raw1 = fastq_scan_raw_reads.read1_seq
    Int fastq_scan_num_reads_clean1 = fastq_scan_clean_reads.read1_seq
    String fastq_scan_version = fastq_scan_raw_reads.version
    File fastq_scan_raw1_json = fastq_scan_raw_reads.fastq_scan_json
    File fastq_scan_clean1_json = fastq_scan_clean_reads.fastq_scan_json
    # Read QC - kraken outputs
    String kraken_version = kraken2_raw.version
    Float kraken_human = kraken2_raw.percent_human
    String kraken_sc2 = kraken2_raw.percent_sc2
    String kraken_target_organism = kraken2_raw.percent_target_organism
    String kraken_target_organism_name = organism_parameters.kraken_target_organism
    File kraken_report = kraken2_raw.kraken_report
    Float kraken_human_dehosted = kraken2_dehosted.percent_human
    String kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String kraken_target_organism_dehosted = kraken2_dehosted.percent_target_organism
    File kraken_report_dehosted = kraken2_dehosted.kraken_report
    # Read Alignment - Artic consensus outputs
    File aligned_bam = consensus.trim_sorted_bam
    File aligned_bai = consensus.trim_sorted_bai
    File variants_from_ref_vcf = consensus.medaka_pass_vcf
    File assembly_fasta = consensus.consensus_seq
    File? read1_aligned = consensus.reads_aligned
    # Read Alignment - Artic consensus versioning outputs
    String artic_version = consensus.artic_pipeline_version
    String artic_docker = consensus.artic_pipeline_docker
    String medaka_reference = consensus.medaka_reference
    String primer_bed_name = consensus.primer_bed_name
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): ~{consensus.artic_pipeline_version}"
    # Read Alignment - consensus assembly qc outputs
    File consensus_stats = stats_n_coverage.stats
    File consensus_flagstat = stats_n_coverage.flagstat
    Float meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float assembly_mean_coverage = stats_n_coverage_primtrim.depth
    String samtools_version_stats = stats_n_coverage.samtools_version
    # Read Alignment - consensus assembly summary outputs
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Percentage mapped reads
    Float percentage_mapped_reads = stats_n_coverage.percentage_mapped_reads
    # SC2 specific coverage outputs
    Float? sc2_s_gene_mean_coverage = gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = gene_coverage.sc2_s_gene_percent_coverage
    File? est_percent_gene_coverage_tsv = gene_coverage.est_percent_gene_coverage_tsv
    # Pangolin outputs
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report= pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs
    File? nextclade_json = nextclade_v3.nextclade_json
    File? auspice_json = nextclade_v3.auspice_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    String? nextclade_version = nextclade_v3.nextclade_version
    String? nextclade_docker = nextclade_v3.nextclade_docker
    String nextclade_ds_tag = organism_parameters.nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
  }
}