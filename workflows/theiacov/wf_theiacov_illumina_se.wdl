version 1.0

import "../utilities/wf_read_QC_trim_se.wdl" as read_qc
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/gene_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_illumina_se {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File? primer_bed
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-02-25T12:00:00Z"
    String? nextclade_dataset_name
    File? reference_genome
    Int min_depth = 100
    String organism = "sars-cov-2"
    Boolean trim_primers = true
    File? adapters
    File? phix
  }
  call read_qc.read_QC_trim_se as read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1_raw,
      adapters = adapters,
      phix = phix,
      workflow_series = "theiacov"
  }
  call consensus_call.ivar_consensus {
    input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      reference_genome = reference_genome,
      primer_bed = primer_bed,
      min_depth = min_depth,
      trim_primers = trim_primers
    }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = ivar_consensus.assembly_fasta,
      reference_genome = reference_genome
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = ivar_consensus.assembly_fasta
    }
    call sc2_calculation.sc2_gene_coverage {
      input: 
        samplename = samplename,
        bamfile = ivar_consensus.aligned_bam,
        min_depth = min_depth
    }
  }
  if (organism == "MPXV") {
    # MPXV specific tasks
  }
  if (organism == "WNV") {
    # WNV specific tasks (none yet, just adding as placeholder for future)
  }
  if (organism == "MPXV" || organism == "sars-cov-2"){
    # tasks specific to either MPXV or sars-cov-2
    call nextclade_task.nextclade {
      input:
      genome_fasta = ivar_consensus.assembly_fasta,
      dataset_name = select_first([nextclade_dataset_name, organism,]),
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
    }
    call nextclade_task.nextclade_output_parser {
      input:
      nextclade_tsv = nextclade.nextclade_tsv,
      organism = organism
    }
  }
  if (organism == "MPXV" || organism == "sars-cov-2" || organism == "WNV"){ 
    # tasks specific to MPXV, sars-cov-2, and WNV
    call vadr_task.vadr {
      input:
        genome_fasta = ivar_consensus.assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_se_version = version_capture.phb_version
    String theiacov_illumina_se_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC - fastq_scan outputs
    Int num_reads_raw = read_QC_trim.fastq_scan_raw_number_reads
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    Int num_reads_clean = read_QC_trim.fastq_scan_clean_number_reads
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    # Read QC - bbduk outputs
    File read1_clean = read_QC_trim.read1_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - kraken outputs
    Float? kraken_human = read_QC_trim.kraken_human
    Float? kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_org = read_QC_trim.kraken_target_org
    String? kraken_target_org_name = read_QC_trim.kraken_target_org_name
    String? kraken_version = read_QC_trim.kraken_version
    File? kraken_report = read_QC_trim.kraken_report
    # Read Alignment - bwa outputs
    String bwa_version = ivar_consensus.bwa_version
    String samtools_version = ivar_consensus.samtools_version
    File read1_aligned = ivar_consensus.read1_aligned
    String assembly_method = ivar_consensus.assembly_method_nonflu
    File aligned_bam = ivar_consensus.aligned_bam
    File aligned_bai = ivar_consensus.aligned_bai
    # Read Alignment - primer trimming outputs
    Float? primer_trimmed_read_percent = ivar_consensus.primer_trimmed_read_percent
    String? ivar_version_primtrim = ivar_consensus.ivar_version_primtrim
    String? samtools_version_primtrim = ivar_consensus.samtools_version_primtrim
    String? primer_bed_name = ivar_consensus.primer_bed_name
    # Read Alignment - variant call outputs
    File ivar_tsv = ivar_consensus.ivar_tsv
    File ivar_vcf = ivar_consensus.ivar_vcf
    String ivar_variant_version = ivar_consensus.ivar_variant_version
    # Read Alignment - assembly outputs
    File assembly_fasta = ivar_consensus.assembly_fasta
    String ivar_version_consensus = ivar_consensus.ivar_version_consensus
    String samtools_version_consensus = ivar_consensus.samtools_version_consensus
    # Read Alignment - consensus assembly qc outputs
    Int consensus_n_variant_min_depth = min_depth
    File consensus_stats = ivar_consensus.consensus_stats
    File consensus_flagstat = ivar_consensus.consensus_flagstat
    Float meanbaseq_trim = ivar_consensus.meanbaseq_trim
    Float meanmapq_trim = ivar_consensus.meanmapq_trim
    Float assembly_mean_coverage = ivar_consensus.assembly_mean_coverage
    String samtools_version_stats = ivar_consensus.samtools_version_stats
    # Read Alignment - consensus assembly summary outputs
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # SC2 specific coverage outputs
    Float? sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage
    File? sc2_all_genes_percent_coverage = sc2_gene_coverage.sc2_all_genes_percent_coverage
    # Pangolin outputs
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs
    File? nextclade_json = nextclade.nextclade_json
    File? auspice_json = nextclade.auspice_json
    File? nextclade_tsv = nextclade.nextclade_tsv
    String? nextclade_version = nextclade.nextclade_version
    String? nextclade_docker = nextclade.nextclade_docker
    String nextclade_ds_tag = nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
  }
}