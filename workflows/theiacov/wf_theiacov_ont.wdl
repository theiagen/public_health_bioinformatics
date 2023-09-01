version 1.0

import "../../tasks/assembly/task_artic_consensus.wdl" as artic_consensus
import "../../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/task_screen.wdl" as screen
import "../../tasks/quality_control/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/task_quasitools.wdl" as quasitools
import "../../tasks/gene_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning
import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc_trim_workflow

workflow theiacov_ont {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data generated on ONT NGS platforms."
  }
  input {
    String samplename
    File demultiplexed_reads
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method = "OXFORD_NANOPORE"
    File primer_bed
    # assembly parameters
    Int normalise = 200
    Int max_length = 700
    Int min_length = 400
    # nextclade inputs
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-08-17T12:00:00Z"
    String? nextclade_dataset_name
    # reference values
    File? reference_genome
    Int? genome_length
    # kraken inputs
    String? target_org
    # read screen parameters
    Int min_reads = 113 # min basepairs / 300 (which is the longest available read length of an Illumina product)
    Int min_basepairs = 34000 # 20x coverage of hepatitis delta virus
    Int min_genome_size = 1700 # size of hepatitis delta virus
    Int max_genome_size = 2673870 # size of Pandoravirus salinus + 200 kb
    Int min_coverage = 10
    Boolean skip_screen = false
    Boolean skip_mash = false
    # qc check parameters
    File? qc_check_table
  }
  if (organism == "HIV") { # set HIV specific artic version
    String run_prefix = "artic_hiv"
  }
  call screen.check_reads_se as raw_check_reads {
    input:
      read1 = demultiplexed_reads,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_size = min_genome_size,
      max_genome_size = max_genome_size,
      min_coverage = min_coverage,
      skip_screen = skip_screen,
      skip_mash = skip_mash,
      workflow_series = "theiacov",
      organism = organism,
      expected_genome_size = genome_length
  }
  if (raw_check_reads.read_screen == "PASS") {
    call read_qc_trim_workflow.read_QC_trim_ont as read_qc_trim {
      input:
        read1 = demultiplexed_reads,
        samplename = samplename,
        genome_size = genome_length,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix,
        target_org = target_org,
        workflow_series = "theiacov"
    }
    call screen.check_reads_se as clean_check_reads {
      input:
        read1 = read_qc_trim.read1_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_size = min_genome_size,
        max_genome_size = max_genome_size,
        min_coverage = min_coverage,
        skip_screen = skip_screen,
        skip_mash = skip_mash,
        workflow_series = "theiacov",
        organism = organism,
        expected_genome_size = genome_length
    }
    if (clean_check_reads.read_screen == "PASS") {
      call artic_consensus.consensus {
        input:
          samplename = samplename,
          organism = organism,
          filtered_reads = read_qc_trim.read1_clean,
          primer_bed = primer_bed,
          normalise = normalise,
          reference_genome = reference_genome
      }
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta = consensus.consensus_seq,
          reference_genome = reference_genome
      }
      # nanoplot for basic QC metrics
      call nanoplot_task.nanoplot as nanoplot_raw {
        input:
          read1 = demultiplexed_reads,
          samplename = samplename,
          est_genome_size = select_first([genome_length, consensus_qc.number_Total])
      }
      call nanoplot_task.nanoplot as nanoplot_clean {
        input:
          read1 = read_qc_trim.read1_clean,
          samplename = samplename,
          est_genome_size = select_first([genome_length, consensus_qc.number_Total])
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
        call nextclade_task.nextclade {
          input:
          genome_fasta = consensus.consensus_seq,
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
            genome_fasta = consensus.consensus_seq,
            assembly_length_unambiguous = consensus_qc.number_ATCG
        }
      }
      if (organism == "HIV") {
        call quasitools.quasitools as quasitools_ont {
          input:
            read1 = read_qc_trim.read1_clean,
            samplename = samplename
        }
      }
      if (defined(qc_check_table)) {
        call qc_check.qc_check_phb as qc_check_task {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = organism,
            num_reads_raw1 = nanoplot_raw.num_reads,
            num_reads_clean1 = nanoplot_clean.num_reads,
            kraken_human = read_qc_trim.kraken_human,
            # kraken_sc2 = kraken2_raw.percent_sc2,
            # kraken_target_org = kraken2_raw.percent_target_org,
            # kraken_human_dehosted = read_QC_trim.kraken_human_dehosted,
            # kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted,
            # kraken_target_org_dehosted =read_QC_trim.kraken_target_org_dehosted,
            meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq,
            assembly_mean_coverage = stats_n_coverage_primtrim.depth,
            number_N = consensus_qc.number_N,
            assembly_length_unambiguous = consensus_qc.number_ATCG,
            number_Degenerate =  consensus_qc.number_Degenerate,
            percent_reference_coverage =  consensus_qc.percent_reference_coverage,
            # sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth,
            # sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage,
            vadr_num_alerts = vadr.num_alerts
        }
      }
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
    # Sample Screening
    String raw_read_screen = raw_check_reads.read_screen
    String? clean_read_screen = clean_check_reads.read_screen
    # Read QC - dehosting outputs
    File? read1_dehosted = read_qc_trim.read1_dehosted
    # Read QC - nanoplot outputs    
    String? nanoplot_version = nanoplot_raw.nanoplot_version
    String? nanoplot_docker = nanoplot_raw.nanoplot_docker
    # Read QC - nanoplot raw outputs
    File? nanoplot_html_raw = nanoplot_raw.nanoplot_html
    File? nanoplot_tsv_raw = nanoplot_raw.nanoplot_tsv
    Int? num_reads_raw1 = nanoplot_raw.num_reads
    Float? r1_mean_readlength_raw = nanoplot_raw.mean_readlength
    Float? r1_mean_q_raw = nanoplot_raw.mean_q
    # Read QC - nanoplot clean outputs
    File? nanoplot_html_clean = nanoplot_clean.nanoplot_html
    File? nanoplot_tsv_clean = nanoplot_clean.nanoplot_tsv
    Int? num_reads_clean1 = nanoplot_clean.num_reads
    Float? r1_mean_readlength_clean = nanoplot_clean.mean_readlength
    Float? r1_mean_q_clean = nanoplot_clean.mean_q
    # Read QC - kraken outputs general
    String? kraken_version = read_qc_trim.kraken_version
    String? kraken_target_org_name = read_qc_trim.kraken_target_org_name
    # Read QC - kraken outputs raw
    Float? kraken_human = read_qc_trim.kraken_human
    Float? kraken_sc2 = read_qc_trim.kraken_sc2
    String? kraken_target_org = read_qc_trim.kraken_target_org
    File? kraken_report = read_qc_trim.kraken_report
    # Read QC - kraken outputs dehosted
    Float? kraken_human_dehosted = read_qc_trim.kraken_human_dehosted
    Float? kraken_sc2_dehosted = read_qc_trim.kraken_sc2_dehosted
    String? kraken_target_org_dehosted = read_qc_trim.kraken_target_org_dehosted
    File? kraken_report_dehosted = read_qc_trim.kraken_report_dehosted
    # Read Alignment - Artic consensus outputs
    File? aligned_bam = consensus.trim_sorted_bam
    File? aligned_bai = consensus.trim_sorted_bai
    File? variants_from_ref_vcf = consensus.medaka_pass_vcf
    File? assembly_fasta = consensus.consensus_seq
    File? read1_aligned = consensus.reads_aligned
    File? read1_trimmed = consensus.trim_fastq
    # Read Alignment - Artic consensus versioning outputs
    String? artic_version = consensus.artic_pipeline_version
    String? artic_docker = consensus.artic_pipeline_docker
    String? medaka_reference = consensus.medaka_reference
    String? primer_bed_name = consensus.primer_bed_name
    String? assembly_method = "TheiaCoV (~{version_capture.phb_version}): ~{consensus.artic_pipeline_version}"
    # Assembly QC - consensus assembly qc outputs
    File? consensus_stats = stats_n_coverage.stats
    File? consensus_flagstat = stats_n_coverage.flagstat
    Float? meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float? meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    Float? assembly_mean_coverage = stats_n_coverage_primtrim.depth
    String? samtools_version = stats_n_coverage.samtools_version
    # Assembly QC - consensus assembly summary outputs
    Int? number_N = consensus_qc.number_N
    Int? assembly_length_unambiguous = consensus_qc.number_ATCG
    Int? number_Degenerate = consensus_qc.number_Degenerate
    Int? number_Total = consensus_qc.number_Total
    Float? percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Assembly QC - nanoplot outputs
    Float? est_coverage_raw = nanoplot_raw.est_coverage
    Float? est_coverage_clean = nanoplot_clean.est_coverage
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
    # HIV outputs
    String? quasitools_version = quasitools_ont.quasitools_version
    String? quasitools_date = quasitools_ont.quasitools_date
    File? quasitools_coverage_file = quasitools_ont.coverage_file
    File? quasitools_dr_report = quasitools_ont.dr_report
    File? quasitools_hydra_vcf = quasitools_ont.hydra_vcf
    File? quasitools_mutations_report = quasitools_ont.mutations_report
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
  }
}