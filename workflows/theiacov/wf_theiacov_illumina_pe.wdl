version 1.0

import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/lentivirus/task_quasitools.wdl" as quasitools
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_flu_track.wdl" as run_flu_track
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc

workflow theiacov_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String organism = "sars-cov-2" # recommended options: "sars-cov-2", "HIV", "WNV", "MPXV", "flu", "rsv-a", "rsv-b"
    File read1
    File read2
    # sequencing values
    String seq_method = "ILLUMINA"
    File? primer_bed
    File? adapters
    File? phix
    # reference values
    File? reference_gff
    File? reference_genome
    File? reference_gene_locations_bed
    Int? genome_length 
    # trimming parameters
    Boolean trim_primers = true
    Int trim_min_length = 75
    Int trim_quality_min_score = 30
    Int trim_window_size = 4
    # assembly parameters
    Int? min_depth # minimum depth to use for consensus and variant calling; default is 100 for non-flu (default value set below in call block for ivar consensus subwf), flu default is 30 for illumina (default set below in flu_track call block)
    Float consensus_min_freq = 0.6 # minimum frequency for a variant to be called as SNP in consensus genome
    Float variant_min_freq = 0.6 # minimum frequency for a variant to be reported in ivar outputs
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_options
    File? vadr_model_file
    Int? vadr_memory
    # read screen parameters
    Int min_reads = 57 # min basepairs / 300 (which is the longest available read length of an Illumina product)
    Int min_basepairs = 17000 # 10x coverage of hepatitis delta virus
    Int min_genome_length = 1700 # size of hepatitis delta virus
    Int max_genome_length = 2673870 # size of Pandoravirus salinus + 200 kb
    Int min_coverage = 10
    Int min_proportion = 40
    Boolean skip_screen = false
    # pangolin parameters
    String? pangolin_docker_image
    # kraken parameters
    String? target_organism
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
      pangolin_docker_image = pangolin_docker_image,
      kraken_target_organism_input = target_organism
  }
  if (! skip_screen) {
    call screen.check_reads as raw_check_reads {
      input:
        read1 = read1,
        read2 = read2,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        workflow_series = "theiacov",
        expected_genome_length = organism_parameters.genome_length
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    call read_qc.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2,
        adapters = adapters,
        phix = phix,
        workflow_series = "theiacov",
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        target_organism = organism_parameters.kraken_target_organism
    }
    if (! skip_screen) {
      call screen.check_reads as clean_check_reads {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          min_reads = min_reads,
          min_basepairs = min_basepairs,
          min_genome_length = min_genome_length,
          max_genome_length = max_genome_length,
          min_coverage = min_coverage,
          min_proportion = min_proportion,
          workflow_series = "theiacov",
          expected_genome_length = organism_parameters.genome_length
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      # assembly via bwa and ivar for non-flu data
      if (organism_parameters.standardized_organism != "flu") {
        call consensus_call.ivar_consensus {
          input:
            samplename = samplename,
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            reference_genome = organism_parameters.reference,
            primer_bed = organism_parameters.primer_bed,
            reference_gff = organism_parameters.reference_gff,
            min_depth = select_first([min_depth, 100]),
            consensus_min_freq = consensus_min_freq,
            variant_min_freq = variant_min_freq,
            trim_primers = trim_primers
        }
      }
      # for flu organisms call flu_track
      if (organism_parameters.standardized_organism == "flu") {
        call run_flu_track.flu_track {
          input:
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            samplename = samplename,
            standardized_organism = organism_parameters.standardized_organism,
            seq_method = seq_method,
            irma_min_consensus_support = select_first([min_depth, 30])
        }
      }
      if (defined(ivar_consensus.assembly_fasta) || defined(flu_track.irma_assembly_fasta)) {
        call consensus_qc_task.consensus_qc {
          input:
            assembly_fasta =  select_first([ivar_consensus.assembly_fasta, flu_track.irma_assembly_fasta]),
            reference_genome = organism_parameters.reference,
            genome_length = organism_parameters.genome_length
        }
        # run organism-specific typing
        if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "measles") { 
          # tasks specific to either MPXV, sars-cov-2, or RSV-A/RSV-B
          call nextclade_task.nextclade_v3 {
            input:
              genome_fasta = select_first([ivar_consensus.assembly_fasta]),
              dataset_name = organism_parameters.nextclade_dataset_name,
              dataset_tag = organism_parameters.nextclade_dataset_tag
          }
          call nextclade_task.nextclade_output_parser {
            input:
              nextclade_tsv = nextclade_v3.nextclade_tsv,
              organism = organism_parameters.standardized_organism
          }
        }
        if (organism_parameters.standardized_organism == "sars-cov-2") {
          # sars-cov-2 specific tasks
          call pangolin.pangolin4 {
            input:
              samplename = samplename,
              fasta = select_first([ivar_consensus.assembly_fasta]),
              docker = organism_parameters.pangolin_docker
          }
        }
        if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || defined(reference_gene_locations_bed)) {
          # tasks specific to either sars-cov-2, MPXV, or any organism with a user-supplied reference gene locations bed file
          call gene_coverage_task.gene_coverage {
            input:
              bamfile = select_first([ivar_consensus.aligned_bam, flu_track.irma_ha_bam, flu_track.irma_na_bam, ""]),
              bedfile = select_first([reference_gene_locations_bed, organism_parameters.gene_locations_bed]),
              samplename = samplename,
              organism = organism_parameters.standardized_organism
          }
        }
        if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles") {
          # tasks specific to MPXV, sars-cov-2, WNV, flu, rsv_a, and rsv_b, measles, mumps, and rubella
          call vadr_task.vadr {
            input:
              genome_fasta = select_first([ivar_consensus.assembly_fasta, flu_track.irma_assembly_fasta_padded]),
              assembly_length_unambiguous = consensus_qc.number_ATCG,
              vadr_opts = organism_parameters.vadr_opts,
              vadr_model_file = organism_parameters.vadr_model_file,
              max_length = organism_parameters.vadr_maxlength,
              skip_length = organism_parameters.vadr_skiplength,
              memory = organism_parameters.vadr_memory
          }
        }
      }
      if (organism_parameters.standardized_organism == "HIV") {
        call quasitools.quasitools as quasitools_illumina_pe {
          input:
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            samplename = samplename
        }
      }
      if (defined(qc_check_table)) {
        # empty strings for kraken outputs throw an error so avoid those outputs for now
        call qc_check.qc_check_phb as qc_check_task {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = organism_parameters.standardized_organism,
            num_reads_raw1 = read_QC_trim.fastq_scan_raw1,
            num_reads_raw2 = read_QC_trim.fastq_scan_raw2,
            num_reads_clean1 = read_QC_trim.fastq_scan_clean1,
            num_reads_clean2 = read_QC_trim.fastq_scan_clean2,
            kraken_human = read_QC_trim.kraken_human,
            kraken_human_dehosted = read_QC_trim.kraken_human_dehosted,
            meanbaseq_trim = ivar_consensus.meanbaseq_trim,
            assembly_mean_coverage = ivar_consensus.assembly_mean_coverage,
            number_N = consensus_qc.number_N,
            assembly_length_unambiguous = consensus_qc.number_ATCG,
            number_Degenerate =  consensus_qc.number_Degenerate,
            percent_reference_coverage =  consensus_qc.percent_reference_coverage,
            vadr_num_alerts = vadr.num_alerts
        }
      }
    }
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_pe_version = version_capture.phb_version
    String theiacov_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC - fastq_scan outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker    
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # Read QC - bbduk outputs
    File? read1_clean = read_QC_trim.read1_clean
    File? read2_clean = read_QC_trim.read2_clean
    String? bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
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
    File? read2_aligned = ivar_consensus.read2_aligned
    String aligned_bam = select_first([ivar_consensus.aligned_bam, ""])
    String aligned_bai = select_first([ivar_consensus.aligned_bai, ""])
    File? read1_unaligned = ivar_consensus.read1_unaligned
    File? read2_unaligned = ivar_consensus.read2_unaligned
    File? sorted_bam_unaligned = ivar_consensus.sorted_bam_unaligned
    File? sorted_bam_unaligned_bai = ivar_consensus.sorted_bam_unaligned_bai
    # Read Alignment - primer trimming outputs
    Float? primer_trimmed_read_percent = ivar_consensus.primer_trimmed_read_percent
    String? ivar_version_primtrim = ivar_consensus.ivar_version_primtrim
    String? samtools_version_primtrim = ivar_consensus.samtools_version
    String? primer_bed_name = ivar_consensus.primer_bed_name
    # Read Alignment - variant call outputs
    File? ivar_tsv = ivar_consensus.ivar_tsv
    File? ivar_vcf = ivar_consensus.ivar_vcf
    String? ivar_variant_proportion_intermediate = ivar_consensus.ivar_variant_proportion_intermediate
    String? ivar_variant_version = ivar_consensus.ivar_variant_version
    String? samtools_version_consensus = ivar_consensus.samtools_version_consensus
    # Read Alignment - assembly outputs
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([ivar_consensus.assembly_method_nonflu, flu_track.irma_version, ""])
    String assembly_fasta = select_first([ivar_consensus.assembly_fasta, flu_track.irma_assembly_fasta, "Assembly could not be generated"])
    String? ivar_version_consensus = ivar_consensus.ivar_version_consensus
    # Read Alignment - consensus assembly qc outputs
    # this is the minimum depth used for consensus and variant calling in EITHER iVar or IRMA
    Int consensus_n_variant_min_depth = select_first([min_depth, flu_track.irma_minimum_consensus_support, 100])
    File? consensus_stats = ivar_consensus.consensus_stats
    File? consensus_flagstat = ivar_consensus.consensus_flagstat
    String meanbaseq_trim = select_first([ivar_consensus.meanbaseq_trim, ""])
    String meanmapq_trim = select_first([ivar_consensus.meanmapq_trim, ""])
    String assembly_mean_coverage = select_first([ivar_consensus.assembly_mean_coverage, flu_track.ha_na_assembly_coverage , ""])
    String? samtools_version_stats = ivar_consensus.samtools_version_stats
    # Read Alignment - consensus assembly summary outputs
    Int? number_N = consensus_qc.number_N
    Int? assembly_length_unambiguous = consensus_qc.number_ATCG
    Int? number_Degenerate =  consensus_qc.number_Degenerate
    Int? number_Total = consensus_qc.number_Total
    Float? percent_reference_coverage =  consensus_qc.percent_reference_coverage
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
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs for all organisms
    String nextclade_version = select_first([nextclade_v3.nextclade_version, flu_track.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade_v3.nextclade_docker, flu_track.nextclade_docker, ""])
    # Nextclade outputs for non-flu
    File? nextclade_json = nextclade_v3.nextclade_json
    File? auspice_json = nextclade_v3.auspice_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    String nextclade_ds_tag = organism_parameters.nextclade_dataset_tag
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # Nextclade outputs for flu H5N1
    File? nextclade_json_flu_h5n1 = flu_track.nextclade_json_flu_h5n1
    File? auspice_json_flu_h5n1 = flu_track.auspice_json_flu_h5n1
    File? nextclade_tsv_flu_h5n1 = flu_track.nextclade_tsv_flu_h5n1
    String? nextclade_aa_subs_flu_h5n1 = flu_track.nextclade_aa_subs_flu_h5n1
    String? nextclade_aa_dels_flu_h5n1 = flu_track.nextclade_aa_dels_flu_h5n1
    String? nextclade_clade_flu_h5n1 = flu_track.nextclade_clade_flu_h5n1
    String? nextclade_qc_flu_h5n1 = flu_track.nextclade_qc_flu_h5n1
    # Nextclade outputs for flu HA
    File? nextclade_json_flu_ha = flu_track.nextclade_json_flu_ha
    File? auspice_json_flu_ha = flu_track.auspice_json_flu_ha
    File? nextclade_tsv_flu_ha = flu_track.nextclade_tsv_flu_ha
    String? nextclade_ds_tag_flu_ha = flu_track.nextclade_ds_tag_flu_ha
    String? nextclade_aa_subs_flu_ha = flu_track.nextclade_aa_subs_flu_ha
    String? nextclade_aa_dels_flu_ha = flu_track.nextclade_aa_dels_flu_ha
    String? nextclade_clade_flu_ha = flu_track.nextclade_clade_flu_ha
    String? nextclade_qc_flu_ha = flu_track.nextclade_qc_flu_ha
    # Nextclade outputs for flu NA
    File? nextclade_json_flu_na = flu_track.nextclade_json_flu_na
    File? auspice_json_flu_na = flu_track.auspice_json_flu_na
    File? nextclade_tsv_flu_na = flu_track.nextclade_tsv_flu_na
    String? nextclade_ds_tag_flu_na = flu_track.nextclade_ds_tag_flu_na
    String? nextclade_aa_subs_flu_na = flu_track.nextclade_aa_subs_flu_na
    String? nextclade_aa_dels_flu_na = flu_track.nextclade_aa_dels_flu_na
    String? nextclade_clade_flu_na = flu_track.nextclade_clade_flu_na
    String? nextclade_qc_flu_na = flu_track.nextclade_qc_flu_na
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # Flu IRMA Outputs
    String? irma_version = flu_track.irma_version
    String? irma_docker = flu_track.irma_docker
    String? irma_type = flu_track.irma_type
    String? irma_subtype = flu_track.irma_subtype
    String? irma_subtype_notes = flu_track.irma_subtype_notes
    File? irma_assembly_fasta_concatenated = flu_track.flu_assembly_fasta_concatenated
    File? irma_ha_segment_fasta = flu_track.flu_ha_segment_fasta
    File? irma_na_segment_fasta = flu_track.flu_na_segment_fasta
    File? irma_pa_segment_fasta = flu_track.flu_pa_segment_fasta
    File? irma_pb1_segment_fasta = flu_track.flu_pb1_segment_fasta
    File? irma_pb2_segment_fasta = flu_track.flu_pb2_segment_fasta
    File? irma_mp_segment_fasta = flu_track.flu_mp_segment_fasta
    File? irma_np_segment_fasta = flu_track.flu_np_segment_fasta
    File? irma_ns_segment_fasta = flu_track.flu_ns_segment_fasta
    # Flu GenoFLU Outputs
    String? genoflu_version = flu_track.genoflu_version
    String? genoflu_genotype = flu_track.genoflu_genotype
    String? genoflu_all_segments = flu_track.genoflu_all_segments
    File? genoflu_output_tsv = flu_track.genoflu_output_tsv
    # Flu Abricate Outputs
    String? abricate_flu_type = flu_track.abricate_flu_type
    String? abricate_flu_subtype =  flu_track.abricate_flu_subtype
    File? abricate_flu_results = flu_track.abricate_flu_results
    String? abricate_flu_database =  flu_track.abricate_flu_database
    String? abricate_flu_version = flu_track.abricate_flu_version
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = flu_track.flu_A_315675_resistance
    String? flu_amantadine_resistance = flu_track.flu_amantadine_resistance
    String? flu_compound_367_resistance = flu_track.flu_compound_367_resistance
    String? flu_favipiravir_resistance = flu_track.flu_favipiravir_resistance
    String? flu_fludase_resistance = flu_track.flu_fludase_resistance
    String? flu_L_742_001_resistance = flu_track.flu_L_742_001_resistance
    String? flu_laninamivir_resistance = flu_track.flu_laninamivir_resistance
    String? flu_peramivir_resistance = flu_track.flu_peramivir_resistance
    String? flu_pimodivir_resistance = flu_track.flu_pimodivir_resistance
    String? flu_rimantadine_resistance = flu_track.flu_rimantadine_resistance
    String? flu_oseltamivir_resistance = flu_track.flu_oseltamivir_resistance
    String? flu_xofluza_resistance = flu_track.flu_xofluza_resistance
    String? flu_zanamivir_resistance = flu_track.flu_zanamivir_resistance
    # HIV Outputs
    String? quasitools_version = quasitools_illumina_pe.quasitools_version
    String? quasitools_date = quasitools_illumina_pe.quasitools_date
    File? quasitools_coverage_file = quasitools_illumina_pe.coverage_file
    File? quasitools_dr_report = quasitools_illumina_pe.dr_report
    File? quasitools_hydra_vcf = quasitools_illumina_pe.hydra_vcf
    File? quasitools_mutations_report = quasitools_illumina_pe.mutations_report
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # Capture percentage_mapped_reads from ivar_consensus task or flu_track task
    String percentage_mapped_reads = select_first([ivar_consensus.percentage_mapped_reads, flu_track.percentage_mapped_reads, ""])
  }
}
