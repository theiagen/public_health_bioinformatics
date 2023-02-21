version 1.0

import "../../tasks/assembly/task_artic_consensus.wdl" as artic_consensus
import "../../tasks/quality_control/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken2
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/task_quasitools.wdl" as quasitools
import "../../tasks/gene_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_ont {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data generated on ONT NGS platforms."
  }
  input {
    String samplename
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    File demultiplexed_reads
    Int normalise = 200
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2022-12-14T12:00:00Z"
    String? nextclade_dataset_name
    File? reference_genome
    Int max_length = 700
    Int min_length = 400
    String organism = "sars-cov-2"
    String? target_org
  }
  call fastq_scan.fastq_scan_se as fastq_scan_raw_reads {
    input:
      read1 = demultiplexed_reads
  }
  call ncbi_scrub.ncbi_scrub_se {
    input:
      samplename = samplename,
      read1 = demultiplexed_reads
  }
  call artic_guppyplex.read_filtering {
    input:
      demultiplexed_reads = ncbi_scrub_se.read1_dehosted,
      samplename = samplename,
      min_length = min_length,
      max_length = max_length
  }
  call fastq_scan.fastq_scan_se as fastq_scan_clean_reads {
    input:
      read1 = read_filtering.filtered_reads
  }
  call kraken2.kraken2_theiacov as kraken2_raw {
    input:
      samplename = samplename,
      read1 = demultiplexed_reads,
      target_org = target_org
  }  
  call kraken2.kraken2_theiacov as kraken2_dehosted {
    input:
      samplename = samplename,
      read1 = ncbi_scrub_se.read1_dehosted,
      target_org = target_org
  }
  call artic_consensus.consensus {
    input:
      samplename = samplename,
      organism = organism,
      filtered_reads = read_filtering.filtered_reads,
      primer_bed = primer_bed,
      normalise = normalise,
      reference_genome = reference_genome
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = consensus.consensus_seq,
      reference_genome = reference_genome
  }
  call assembly_metrics.stats_n_coverage {
    input:
      samplename = samplename,
      bamfile = consensus.sorted_bam
  }
  call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
    input:
      samplename = samplename,
      bamfile = consensus.trim_sorted_bam
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = consensus.consensus_seq
    }
    call sc2_calculation.sc2_gene_coverage {
      input: 
        samplename = samplename,
        bamfile = consensus.trim_sorted_bam,
        min_depth = 20
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
    call nextclade.nextclade_one_sample {
      input:
      genome_fasta = consensus.consensus_seq,
      dataset_name = select_first([nextclade_dataset_name, organism,]),
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
    }
    call nextclade.nextclade_output_parser_one_sample {
      input:
      nextclade_tsv = nextclade_one_sample.nextclade_tsv
    }
   }
  if (organism == "MPXV" || organism == "sars-cov-2" || organism == "WNV"){ 
    # tasks specific to MPXV, sars-cov-2, and WNV
    call vadr_task.vadr {
      input:
        genome_fasta = consensus.consensus_seq,
        assembly_length_unambiguous = consensus_qc.number_ATCG
    }
  }
  if (organism == "HIV") {
    call quasitools.quasitools as quasitools_ont {
      input:
        read1 = read_filtering.filtered_reads,
        samplename = samplename
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_ont_version = version_capture.phb_version
    String theiacov_ont_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC
    File read1_dehosted = ncbi_scrub_se.read1_dehosted
    Int num_reads_raw1 = fastq_scan_raw_reads.read1_seq
    Int num_reads_clean1 = fastq_scan_clean_reads.read1_seq
    String fastq_scan_version = fastq_scan_clean_reads.version
    String kraken_version = kraken2_raw.version
    Float kraken_human = kraken2_raw.percent_human
    Float kraken_sc2 = kraken2_raw.percent_sc2
    String? kraken_target_org = kraken2_raw.percent_target_org
    String? kraken_target_org_name = kraken2_raw.kraken_target_org
    File kraken_report = kraken2_raw.kraken_report
    Float kraken_human_dehosted = kraken2_dehosted.percent_human
    Float kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String? kraken_target_org_dehosted = kraken2_dehosted.percent_target_org
    File kraken_report_dehosted = kraken2_dehosted.kraken_report
    # Read Alignment
    File aligned_bam = consensus.trim_sorted_bam
    File aligned_bai = consensus.trim_sorted_bai
    File variants_from_ref_vcf = consensus.medaka_pass_vcf
    String artic_version = consensus.artic_pipeline_version
    String artic_docker = consensus.artic_pipeline_docker
    String medaka_reference = consensus.medaka_reference
    String primer_bed_name = consensus.primer_bed_name
    File assembly_fasta = consensus.consensus_seq
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): ~{consensus.artic_pipeline_version}"

    File? read1_aligned = consensus.reads_aligned
    File? trim_fastq = consensus.trim_fastq
    # Assembly QC
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Alignment QC
    File consensus_stats = stats_n_coverage.stats
    File consensus_flagstat = stats_n_coverage.flagstat
    Float meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float assembly_mean_coverage = stats_n_coverage_primtrim.depth
    String samtools_version = stats_n_coverage.samtools_version
    # SC2 Specific
    Float? sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage
    File? sc2_all_genes_percent_coverage = sc2_gene_coverage.sc2_all_genes_percent_coverage
    # Lineage Assignment
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Clade Assigment
    File? nextclade_json = nextclade_one_sample.nextclade_json
    File? auspice_json = nextclade_one_sample.auspice_json
    File? nextclade_tsv = nextclade_one_sample.nextclade_tsv
    String? nextclade_version = nextclade_one_sample.nextclade_version
    String? nextclade_docker = nextclade_one_sample.nextclade_docker
    String nextclade_ds_tag = nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser_one_sample.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser_one_sample.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser_one_sample.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser_one_sample.nextclade_lineage
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # HIV outputs
    String? quasitools_version = quasitools_ont.quasitools_version
    String? quasitools_date = quasitools_ont.quasitools_date
    File? quasitools_coverage_file = quasitools_ont.coverage_file
    File? quasitools_dr_report = quasitools_ont.dr_report
    File? quasitools_hydra_vcf = quasitools_ont.hydra_vcf
    File? quasitools_mutations_report = quasitools_ont.mutations_report
  }
}