version 1.0

import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/task_versioning.wdl" as versioning
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_read_QC_trim_se.wdl" as read_qc
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf

workflow theiacov_illumina_se {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    File read1
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method = "ILLUMINA"
    File? primer_bed
    File? adapters
    File? phix
    # trimming parameters
    Boolean trim_primers = true
    Int trim_min_length = 25
    Int trim_quality_min_score = 30
    Int trim_window_size = 4
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # reference values
    File? reference_gff
    File? reference_genome
    File? reference_gene_locations_bed
    Int? genome_length
    # assembly parameters
    Int min_depth = 100
    Float consensus_min_freq = 0.6 # minimum frequency for a variant to be called as SNP in consensus genome
    Float variant_min_freq = 0.6 # minimum frequency for a variant to be reported in ivar outputs
    # read screen parameters
    Int min_reads = 57 # min basepairs / 300 (which is the longest available read length of an Illumina product)
    Int min_basepairs = 17000 # 10x coverage of hepatitis delta virus
    Int min_genome_length = 1700 # size of hepatitis delta virus
    Int max_genome_length = 2673870 # size of Pandoravirus salinus + 200 kb
    Int min_coverage = 10
    Boolean skip_screen = false
    Boolean skip_mash = false
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_options
    File? vadr_model_file
    Int? vadr_memory
    # pangolin parameters
    String? pangolin_docker_image
    # qc check parameters
    File? qc_check_table
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      reference_gff_file = reference_gff,
      reference_genome = reference_genome,
      gene_locations_bed_file = reference_gene_locations_bed,
      genome_length_input = genome_length,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,     
      vadr_max_length = vadr_max_length,
      vadr_skip_length = vadr_skip_length,
      vadr_options = vadr_options,
      vadr_model = vadr_model_file,
      vadr_mem = vadr_memory,
      primer_bed_file = primer_bed,
      pangolin_docker_image = pangolin_docker_image  
  }
  if (! skip_screen) {
    call screen.check_reads_se as raw_check_reads {
      input:
        read1 = read1,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        workflow_series = "theiacov",
        skip_mash = skip_mash,
        expected_genome_length = organism_parameters.genome_length
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    call read_qc.read_QC_trim_se as read_QC_trim {
      input:
        samplename = samplename,
        read1 = read1,
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        adapters = adapters,
        phix = phix,
        workflow_series = "theiacov",
        target_organism = organism_parameters.kraken_target_organism
    }
    if (! skip_screen) {
      call screen.check_reads_se as clean_check_reads {
        input:
          read1 = read_QC_trim.read1_clean,
          min_reads = min_reads,
          min_basepairs = min_basepairs,
          min_genome_length = min_genome_length,
          max_genome_length = max_genome_length,
          min_coverage = min_coverage,
          workflow_series = "theiacov",
          skip_mash = skip_mash,
          expected_genome_length = organism_parameters.genome_length
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      call consensus_call.ivar_consensus {
        input:
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          reference_genome = organism_parameters.reference,
          primer_bed = organism_parameters.primer_bed,
          reference_gff = organism_parameters.reference_gff,
          min_depth = min_depth,
          consensus_min_freq = consensus_min_freq,
          variant_min_freq = variant_min_freq,
          trim_primers = trim_primers
      }
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta = ivar_consensus.assembly_fasta,
          reference_genome = organism_parameters.reference
      }
      call morgana_magic_wf.morgana_magic {
        input:
          samplename = samplename,
          assembly_fasta = ivar_consensus.assembly_fasta,
          taxon_name = organism_parameters.standardized_organism,
          read1 = read_QC_trim.read1_clean,
          number_ATCG = consensus_qc.number_ATCG,
          vadr_max_length = organism_parameters.vadr_maxlength,
          vadr_skip_length = organism_parameters.vadr_skiplength,
          vadr_options = organism_parameters.vadr_opts,
          vadr_model_file = organism_parameters.vadr_model_file,
          vadr_memory = organism_parameters.vadr_memory,
          reference_gene_locations_bed = organism_parameters.gene_locations_bed,
          gene_coverage_bam = ivar_consensus.aligned_bam,
          nextclade_dataset_name = organism_parameters.nextclade_dataset_name,
          nextclade_dataset_tag = organism_parameters.nextclade_dataset_tag,
          pangolin_docker_image = organism_parameters.pangolin_docker,
          # Setting flu_track related inputs to default values as they are not utilized in TheiaCov, decreasing external input bloat
          seq_method = "",
          assembly_metrics_cpu = 0,
          assembly_metrics_disk_size = 0,
          assembly_metrics_docker = "",
          assembly_metrics_memory = 0,
          irma_cpu = 0,
          irma_disk_size = 0,
          irma_docker_image = "",        
          irma_keep_ref_deletions = false,
          irma_memory = 0,
          irma_return_aligned_reads = false,
          genoflu_cpu = 0,
          genoflu_disk_size = 0,
          genoflu_docker = "",
          genoflu_memory = 0,
          abricate_flu_cpu = 0,
          abricate_flu_disk_size = 0,
          abricate_flu_docker = "",
          abricate_flu_memory = 0,
          abricate_flu_min_percent_coverage = 0,
          abricate_flu_min_percent_identity = 0,
          flu_track_antiviral_aa_subs = "",
          workflow_type = "theiacov_se"
      }
      if (defined(qc_check_table)) {
        call qc_check.qc_check_phb as qc_check_task {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = organism_parameters.standardized_organism,
            num_reads_raw1 = select_first([read_QC_trim.fastq_scan_raw1, read_QC_trim.fastqc_raw1]),
            num_reads_clean1 = select_first([read_QC_trim.fastq_scan_clean1, read_QC_trim.fastqc_clean1]),
            kraken_human = read_QC_trim.kraken_human,
            meanbaseq_trim = ivar_consensus.meanbaseq_trim,
            assembly_mean_coverage = ivar_consensus.assembly_mean_coverage,
            number_N = consensus_qc.number_N,
            assembly_length_unambiguous = consensus_qc.number_ATCG,
            number_Degenerate =  consensus_qc.number_Degenerate,
            percent_reference_coverage =  consensus_qc.percent_reference_coverage,
            vadr_num_alerts = morgana_magic.vadr_num_alerts
        }
      }
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_se_version = version_capture.phb_version
    String theiacov_illumina_se_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC - fastq_scan raw outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Float? fastq_scan_r1_mean_readlength_raw = read_QC_trim.fastq_scan_r1_mean_readlength_raw
    Float? fastq_scan_r1_mean_q_raw = read_QC_trim.fastq_scan_r1_mean_q_raw
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? fastq_scan_docker = read_QC_trim.fastq_scan_docker
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    # Read QC - fastq_scan clean outputs
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Float? fastq_scan_r1_mean_readlength_clean = read_QC_trim.fastq_scan_r1_mean_readlength_clean
    Float? fastq_scan_r1_mean_q_clean = read_QC_trim.fastq_scan_r1_mean_q_clean
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # Read QC - bbduk outputs
    File? read1_clean = read_QC_trim.read1_clean
    String? bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - kraken outputs
    String? kraken_version = read_QC_trim.kraken_version
    Float? kraken_human = read_QC_trim.kraken_human
    String? kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_organism = read_QC_trim.kraken_target_organism
    String? kraken_target_organism_name = read_QC_trim.kraken_target_organism_name
    File? kraken_report = read_QC_trim.kraken_report
    Float? kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    String? kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String? kraken_target_organism_dehosted = read_QC_trim.kraken_target_organism_dehosted
    File? kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment - bwa outputs
    String? bwa_version = ivar_consensus.bwa_version
    String? samtools_version = ivar_consensus.samtools_version
    File? read1_aligned = ivar_consensus.read1_aligned
    String? assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([ivar_consensus.assembly_method_nonflu, ""])
    File? aligned_bam = ivar_consensus.aligned_bam
    File? aligned_bai = ivar_consensus.aligned_bai
    File? read1_unaligned = ivar_consensus.read1_unaligned
    File? sorted_bam_unaligned = ivar_consensus.sorted_bam_unaligned
    File? sorted_bam_unaligned_bai = ivar_consensus.sorted_bam_unaligned_bai
    # Read Alignment - primer trimming outputs
    Float? primer_trimmed_read_percent = ivar_consensus.primer_trimmed_read_percent
    String? ivar_version_primtrim = ivar_consensus.ivar_version_primtrim
    String? samtools_version_primtrim = ivar_consensus.samtools_version_primtrim
    String? primer_bed_name = ivar_consensus.primer_bed_name
    # Read Alignment - variant call outputs
    File? ivar_tsv = ivar_consensus.ivar_tsv
    File? ivar_vcf = ivar_consensus.ivar_vcf
    String? ivar_variant_proportion_intermediate = ivar_consensus.ivar_variant_proportion_intermediate
    String? ivar_variant_version = ivar_consensus.ivar_variant_version
    String? samtools_version_consensus = ivar_consensus.samtools_version_consensus
    # Read Alignment - assembly outputs
    File? assembly_fasta = ivar_consensus.assembly_fasta
    String? ivar_version_consensus = ivar_consensus.ivar_version_consensus
    # Read Alignment - consensus assembly qc outputs
    Int? consensus_n_variant_min_depth = min_depth
    File? consensus_stats = ivar_consensus.consensus_stats
    File? consensus_flagstat = ivar_consensus.consensus_flagstat
    Float? meanbaseq_trim = ivar_consensus.meanbaseq_trim
    Float? meanmapq_trim = ivar_consensus.meanmapq_trim
    Float? assembly_mean_coverage = ivar_consensus.assembly_mean_coverage
    String? samtools_version_stats = ivar_consensus.samtools_version_stats
    # Read Alignment - consensus assembly summary outputs
    Int? number_N = consensus_qc.number_N
    Int? assembly_length_unambiguous = consensus_qc.number_ATCG
    Int? number_Degenerate = consensus_qc.number_Degenerate
    Int? number_Total = consensus_qc.number_Total
    Float? percent_reference_coverage = consensus_qc.percent_reference_coverage
    # SC2 specific coverage outputs
    Float? sc2_s_gene_mean_coverage = morgana_magic.sc2_s_gene_mean_coverage
    Float? sc2_s_gene_percent_coverage = morgana_magic.sc2_s_gene_percent_coverage
    File? est_percent_gene_coverage_tsv = morgana_magic.est_percent_gene_coverage_tsv
    # Pangolin outputs
    String? pango_lineage = morgana_magic.pango_lineage
    String? pango_lineage_expanded = morgana_magic.pango_lineage_expanded
    String? pangolin_conflicts = morgana_magic.pangolin_conflicts
    String? pangolin_notes = morgana_magic.pangolin_notes
    String? pangolin_assignment_version = morgana_magic.pangolin_assignment_version
    File? pango_lineage_report = morgana_magic.pango_lineage_report
    String? pangolin_docker = morgana_magic.pangolin_docker
    String? pangolin_versions = morgana_magic.pangolin_versions
    # Nextclade outputs
    File? nextclade_json = morgana_magic.nextclade_json
    File? auspice_json = morgana_magic.auspice_json
    File? nextclade_tsv = morgana_magic.nextclade_tsv
    String? nextclade_version = morgana_magic.nextclade_version
    String? nextclade_docker = morgana_magic.nextclade_docker
    String? nextclade_ds_tag = organism_parameters.nextclade_dataset_tag
    String? nextclade_aa_subs = morgana_magic.nextclade_aa_subs
    String? nextclade_aa_dels = morgana_magic.nextclade_aa_dels
    String? nextclade_clade = morgana_magic.nextclade_clade
    String? nextclade_lineage = morgana_magic.nextclade_lineage
    String? nextclade_qc = morgana_magic.nextclade_qc
    # VADR Annotation QC
    File? vadr_alerts_list = morgana_magic.vadr_alerts_list
    File? vadr_feature_tbl_pass = morgana_magic.vadr_feature_tbl_pass
    File? vadr_feature_tbl_fail = morgana_magic.vadr_feature_tbl_fail
    File? vadr_classification_summary_file = morgana_magic.vadr_classification_summary_file
    File? vadr_all_outputs_tar_gz = morgana_magic.vadr_all_outputs_tar_gz
    String? vadr_num_alerts = morgana_magic.vadr_num_alerts
    String? vadr_docker = morgana_magic.vadr_docker
    File? vadr_fastas_zip_archive = morgana_magic.vadr_fastas_zip_archive
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # Capture percentage_mapped_reads from ivar_consensus task
    String? percentage_mapped_reads = ivar_consensus.percentage_mapped_reads
  }
}