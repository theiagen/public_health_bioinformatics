version 1.0

import "../../tasks/assembly/task_artic_consensus.wdl" as artic_consensus
import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/lentivirus/task_quasitools.wdl" as quasitools
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_flu_track.wdl" as run_flu_track
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc_trim_workflow

workflow theiacov_ont {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data generated on ONT NGS platforms."
  }
  input {
    String samplename
    File read1
    String organism = "sars-cov-2" # options: "sars-cov-2", "HIV", "flu"
    # sequencing values
    String seq_method = "OXFORD_NANOPORE"
    File? primer_bed
    # assembly parameters - sars-cov-2 specific
    Int normalise = 200
    Int max_length = 700
    Int min_length = 400
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # reference values
    File? reference_genome
    File? reference_gene_locations_bed
    Int? genome_length
    # kraken inputs
    String? target_organism
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
    ## flu specific inputs
    # default set to 50 for ONT data in call block below, following CDC MIRA standards
    Int? irma_min_consensus_support
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      reference_genome = reference_genome,
      gene_locations_bed_file = reference_gene_locations_bed,
      genome_length_input = genome_length,
      kraken_target_organism_input = target_organism,
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
  if (organism_parameters.standardized_organism == "HIV") { # set HIV specific artic version
    String run_prefix = "artic_hiv"
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
        skip_mash = skip_mash,
        workflow_series = "theiacov",
        expected_genome_length = organism_parameters.genome_length
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    call read_qc_trim_workflow.read_QC_trim_ont as read_QC_trim {
      input:
        read1 = read1,
        samplename = samplename,
        genome_length = organism_parameters.genome_length,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix,
        target_organism = organism_parameters.kraken_target_organism,
        workflow_series = "theiacov"
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
          skip_mash = skip_mash,
          workflow_series = "theiacov",
            expected_genome_length = organism_parameters.genome_length
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      # assembly via artic_consensus for sars-cov-2 and HIV
      if (organism_parameters.standardized_organism != "flu") {
        call artic_consensus.consensus {
          input:
            samplename = samplename,
            organism = organism_parameters.standardized_organism,
            read1 = read_QC_trim.read1_clean,
            primer_bed = organism_parameters.primer_bed,
            normalise = normalise,
            reference_genome = organism_parameters.reference,
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
      }
      # assembly via irma for flu organisms
      if (organism_parameters.standardized_organism == "flu") {
        call run_flu_track.flu_track {
          input:
            read1 = read_QC_trim.read1_clean,
            samplename = samplename,
            standardized_organism = organism_parameters.standardized_organism,
            seq_method = seq_method,
            irma_min_consensus_support = select_first([irma_min_consensus_support, 50])
        }
      }
      # nanoplot for basic QC metrics
      call nanoplot_task.nanoplot as nanoplot_raw {
        input:
          read1 = read1,
          samplename = samplename,
          est_genome_length = select_first([genome_length, organism_parameters.genome_length])
      }
      call nanoplot_task.nanoplot as nanoplot_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename,
          est_genome_length = select_first([genome_length, organism_parameters.genome_length])
      }
      # consensus QC check
      if (defined(flu_track.irma_assembly_fasta) || defined(consensus.consensus_seq)) {
        call consensus_qc_task.consensus_qc {
          input:
            assembly_fasta =  select_first([flu_track.irma_assembly_fasta, consensus.consensus_seq]),
            reference_genome = organism_parameters.reference,
            genome_length = organism_parameters.genome_length
        }      
        # run organism-specific typing
        if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b") { 
          # tasks specific to either MPXV, sars-cov-2, rsv_a, or rsv_b
          call nextclade_task.nextclade_v3 {
            input:
            genome_fasta = select_first([consensus.consensus_seq]),
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
              fasta = select_first([consensus.consensus_seq]),
              docker = organism_parameters.pangolin_docker
          }
        }  
        if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || defined(reference_gene_locations_bed)) {
          # tasks specific to either sars-cov-2, MPXV, or any organism with a user-supplied reference gene locations bed file
          call gene_coverage_task.gene_coverage {
            input:
              bamfile = select_first([consensus.trim_sorted_bam, flu_track.irma_ha_bam, flu_track.irma_na_bam, ""]),
              bedfile = select_first([reference_gene_locations_bed, organism_parameters.gene_locations_bed]),
              samplename = samplename,
              organism = organism_parameters.standardized_organism
          }
        }
        if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles") {
          # tasks specific to MPXV, sars-cov-2, WNV, flu, rsv_a, and rsv_b, measles, mumps, and rubella
          call vadr_task.vadr {
            input:
              genome_fasta = select_first([consensus.consensus_seq, flu_track.irma_assembly_fasta_padded]),
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
        call quasitools.quasitools as quasitools_ont {
          input:
            read1 = read_QC_trim.read1_clean,
            samplename = samplename
        }
      }
      if (defined(qc_check_table)) {
        call qc_check.qc_check_phb as qc_check_task {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = organism_parameters.standardized_organism,
            num_reads_raw1 = nanoplot_raw.num_reads,
            num_reads_clean1 = nanoplot_clean.num_reads,
            kraken_human = read_QC_trim.kraken_human,
            meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq,
            assembly_mean_coverage = stats_n_coverage_primtrim.depth,
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
    String theiacov_ont_version = version_capture.phb_version
    String theiacov_ont_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    # Read QC - nanoplot outputs    
    String? nanoplot_version = nanoplot_raw.nanoplot_version
    String? nanoplot_docker = nanoplot_raw.nanoplot_docker
    # Read QC - nanoplot raw outputs
    File? nanoplot_html_raw = nanoplot_raw.nanoplot_html
    File? nanoplot_tsv_raw = nanoplot_raw.nanoplot_tsv
    Int? nanoplot_num_reads_raw1 = nanoplot_raw.num_reads
    Float? nanoplot_r1_median_readlength_raw = nanoplot_raw.median_readlength
    Float? nanoplot_r1_mean_readlength_raw = nanoplot_raw.mean_readlength
    Float? nanoplot_r1_stdev_readlength_raw = nanoplot_raw.stdev_readlength
    Float? nanoplot_r1_n50_raw = nanoplot_raw.n50
    Float? nanoplot_r1_mean_q_raw = nanoplot_raw.mean_q
    Float? nanoplot_r1_median_q_raw = nanoplot_raw.median_q
    Float? nanoplot_r1_est_coverage_raw = nanoplot_raw.est_coverage
    # Read QC - nanoplot clean outputs
    File? nanoplot_html_clean = nanoplot_clean.nanoplot_html
    File? nanoplot_tsv_clean = nanoplot_clean.nanoplot_tsv
    Int? nanoplot_num_reads_clean1 = nanoplot_clean.num_reads
    Float? nanoplot_r1_median_readlength_clean = nanoplot_clean.median_readlength
    Float? nanoplot_r1_mean_readlength_clean = nanoplot_clean.mean_readlength
    Float? nanoplot_r1_stdev_readlength_clean = nanoplot_clean.stdev_readlength
    Float? nanoplot_r1_n50_clean = nanoplot_clean.n50
    Float? nanoplot_r1_mean_q_clean = nanoplot_clean.mean_q
    Float? nanoplot_r1_median_q_clean = nanoplot_clean.median_q
    Float? nanoplot_r1_est_coverage_clean = nanoplot_clean.est_coverage
    # Read QC - kraken outputs general
    String? kraken_version = read_QC_trim.kraken_version
    String? kraken_target_organism_name = read_QC_trim.kraken_target_organism_name
    # Read QC - kraken outputs raw
    Float? kraken_human = read_QC_trim.kraken_human
    String? kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_organism = read_QC_trim.kraken_target_organism
    File? kraken_report = read_QC_trim.kraken_report
    # Read QC - kraken outputs dehosted
    Float? kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    String? kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String? kraken_target_organism_dehosted = read_QC_trim.kraken_target_organism_dehosted
    File? kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment - Artic consensus outputs
    String assembly_fasta = select_first([consensus.consensus_seq, flu_track.irma_assembly_fasta, "Assembly could not be generated"])
    File? aligned_bam = consensus.trim_sorted_bam
    File? aligned_bai = consensus.trim_sorted_bai
    File? medaka_vcf = consensus.medaka_pass_vcf
    File? read1_aligned = consensus.reads_aligned
    File? read1_trimmed = consensus.trim_fastq
    # Read Alignment - Artic consensus versioning outputs
    String? artic_version = consensus.artic_pipeline_version
    String? artic_docker = consensus.artic_pipeline_docker
    String? medaka_reference = consensus.medaka_reference
    String? primer_bed_name = consensus.primer_bed_name
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([consensus.artic_pipeline_version, flu_track.irma_version, ""])
    # Assembly QC - consensus assembly qc outputs
    File? consensus_stats = stats_n_coverage.stats
    File? consensus_flagstat = stats_n_coverage.flagstat
    Float? meanbaseq_trim = stats_n_coverage_primtrim.meanbaseq
    Float? meanmapq_trim = stats_n_coverage_primtrim.meanmapq
    String assembly_mean_coverage = select_first([stats_n_coverage_primtrim.depth, flu_track.ha_na_assembly_coverage, ""])
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
    String? vadr_num_alerts = vadr.num_alerts
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # Flu IRMA Outputs
    String? irma_version = flu_track.irma_version
    String? irma_docker = flu_track.irma_docker
    Int? irma_min_consensus_support_threshold = flu_track.irma_minimum_consensus_support
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
    # Non-flu specific outputs
    String percentage_mapped_reads = select_first([stats_n_coverage_primtrim.percentage_mapped_reads, stats_n_coverage.percentage_mapped_reads, flu_track.percentage_mapped_reads, ""])
  }
}