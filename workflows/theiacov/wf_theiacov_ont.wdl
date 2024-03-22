version 1.0

import "../../tasks/assembly/task_artic_consensus.wdl" as artic_consensus
import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate
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
import "../../workflows/utilities/wf_influenza_antiviral_substitutions.wdl" as flu_antiviral
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
    # assembly parameters
    Int normalise = 200
    Int max_length = 700
    Int min_length = 400
    Int min_depth = 20
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
    String? vadr_options
    # pangolin parameters
    String? pangolin_docker_image
    # qc check parameters
    File? qc_check_table
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      reference_genome = reference_genome,
      gene_locations_bed_file = reference_gene_locations_bed,
      genome_length_input = genome_length,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,     
      vadr_max_length = vadr_max_length,
      vadr_options = vadr_options,
      primer_bed_file = primer_bed,
      pangolin_docker_image = pangolin_docker_image
  }
  if (organism_parameters.standardized_organism == "HIV") { # set HIV specific artic version
    String run_prefix = "artic_hiv"
  }
  call screen.check_reads_se as raw_check_reads {
    input:
      read1 = read1,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_length = min_genome_length,
      max_genome_length = max_genome_length,
      min_coverage = min_coverage,
      skip_screen = skip_screen,
      skip_mash = skip_mash,
      workflow_series = "theiacov",
      organism = organism_parameters.standardized_organism,
      expected_genome_length = genome_length
  }
  if (raw_check_reads.read_screen == "PASS") {
    call read_qc_trim_workflow.read_QC_trim_ont as read_qc_trim {
      input:
        read1 = read1,
        samplename = samplename,
        genome_length = genome_length,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix,
        target_organism = organism_parameters.kraken_target_organism,
        workflow_series = "theiacov"
    }
    call screen.check_reads_se as clean_check_reads {
      input:
        read1 = read_qc_trim.read1_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        skip_screen = skip_screen,
        skip_mash = skip_mash,
        workflow_series = "theiacov",
        organism = organism_parameters.standardized_organism,
        expected_genome_length = genome_length
    }
    if (clean_check_reads.read_screen == "PASS") {
      # assembly via artic_consensus for sars-cov-2 and HIV
      if (organism_parameters.standardized_organism != "flu") {
        call artic_consensus.consensus {
          input:
            samplename = samplename,
            organism = organism_parameters.standardized_organism,
            read1 = read_qc_trim.read1_clean,
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
        call irma_task.irma {
          input:
            read1 = read_qc_trim.read1_clean,
            samplename = samplename,
            seq_method = seq_method
        }
        if (defined(irma.irma_assemblies)) {
          call abricate.abricate_flu {
            input:
              assembly = select_first([irma.irma_assembly_fasta]),
              samplename = samplename
          } 
          call set_organism_defaults.organism_parameters as set_flu_na_nextclade_values {
            input:
              organism = organism,
              flu_segment = "NA",
              flu_subtype = irma.irma_subtype,
              # including these to block from terra
              reference_genome = reference_genome,
              genome_length_input = genome_length,
              nextclade_dataset_tag_input = nextclade_dataset_tag,
              nextclade_dataset_name_input = nextclade_dataset_name,     
              vadr_max_length = vadr_max_length,
              vadr_options = vadr_options,
              primer_bed_file = primer_bed,
              pangolin_docker_image = pangolin_docker_image,
              kraken_target_organism_input = target_organism,
              hiv_primer_version = "N/A"
          }
          call set_organism_defaults.organism_parameters as set_flu_ha_nextclade_values {
            input:
              organism = organism,
              flu_segment = "HA",
              flu_subtype = irma.irma_subtype,
               # including these to block from terra
              reference_genome = reference_genome,
              genome_length_input = genome_length,
              nextclade_dataset_tag_input = nextclade_dataset_tag,
              nextclade_dataset_name_input = nextclade_dataset_name,     
              vadr_max_length = vadr_max_length,
              vadr_options = vadr_options,
              primer_bed_file = primer_bed,
              pangolin_docker_image = pangolin_docker_image,
              kraken_target_organism_input = target_organism,
              hiv_primer_version = "N/A"
          }         
          if (set_flu_na_nextclade_values.nextclade_dataset_tag == "NA") {
            Boolean do_not_run_flu_na_nextclade = true
          }
          if (set_flu_ha_nextclade_values.nextclade_dataset_tag == "NA") {
            Boolean do_not_run_flu_ha_nextclade = true
          }
        }
      }
      # consensus QC check
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta =  select_first([irma.irma_assembly_fasta, consensus.consensus_seq]),
          reference_genome = organism_parameters.reference,
          genome_length = organism_parameters.genome_length
      }
      # nanoplot for basic QC metrics
      call nanoplot_task.nanoplot as nanoplot_raw {
        input:
          read1 = read1,
          samplename = samplename,
          est_genome_length = select_first([genome_length, consensus_qc.number_Total, organism_parameters.genome_length])
      }
      call nanoplot_task.nanoplot as nanoplot_clean {
        input:
          read1 = read_qc_trim.read1_clean,
          samplename = samplename,
          est_genome_length = select_first([genome_length, consensus_qc.number_Total, organism_parameters.genome_length])
      }
      if (organism_parameters.standardized_organism == "flu") {
        call flu_antiviral.flu_antiviral_substitutions {
          input:
            na_segment_assembly = irma.seg_na_assembly,
            ha_segment_assembly = irma.seg_ha_assembly,
            pa_segment_assembly = irma.seg_pa_assembly,
            pb1_segment_assembly = irma.seg_pb1_assembly,
            pb2_segment_assembly = irma.seg_pb2_assembly,
            mp_segment_assembly = irma.seg_mp_assembly,
            abricate_flu_subtype = select_first([abricate_flu.abricate_flu_subtype, ""]),
            irma_flu_subtype = select_first([irma.irma_subtype, ""]),
        }
      }
      # run organism-specific typing
      if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || (organism_parameters.standardized_organism == "flu" && defined(irma.seg_ha_assembly) && ! defined(do_not_run_flu_ha_nextclade))) { 
        # tasks specific to either MPXV, sars-cov-2, or flu
        call nextclade_task.nextclade_v3 {
          input:
          genome_fasta = select_first([irma.seg_ha_assembly, consensus.consensus_seq]),
          dataset_name = select_first([set_flu_ha_nextclade_values.nextclade_dataset_name, organism_parameters.nextclade_dataset_name]),
          dataset_tag = select_first([set_flu_ha_nextclade_values.nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag])
        }
        call nextclade_task.nextclade_output_parser {
          input:
          nextclade_tsv = nextclade_v3.nextclade_tsv,
          organism = organism
        }
      }
      if (organism_parameters.standardized_organism == "flu" && defined(irma.seg_na_assembly) && ! defined(do_not_run_flu_na_nextclade)) { 
        # tasks specific to flu NA - run nextclade a second time
        call nextclade_task.nextclade_v3 as nextclade_flu_na {
          input:
            genome_fasta = select_first([irma.seg_na_assembly]),
            dataset_name = select_first([set_flu_na_nextclade_values.nextclade_dataset_name, organism_parameters.nextclade_dataset_name]),
            dataset_tag = select_first([set_flu_na_nextclade_values.nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag])
        }
        call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_na {
          input:
            nextclade_tsv = nextclade_flu_na.nextclade_tsv,
            organism = organism,
            NA_segment = true
        }
        # concatenate tag, aa subs and aa dels for HA and NA segments
        String ha_na_nextclade_ds_tag= "~{set_flu_ha_nextclade_values.nextclade_dataset_tag + ',' + set_flu_na_nextclade_values.nextclade_dataset_tag}"
        String ha_na_nextclade_aa_subs= "~{nextclade_output_parser.nextclade_aa_subs + ',' + nextclade_output_parser_flu_na.nextclade_aa_subs}"
        String ha_na_nextclade_aa_dels= "~{nextclade_output_parser.nextclade_aa_dels + ',' + nextclade_output_parser_flu_na.nextclade_aa_dels}"
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
            bamfile = select_first([consensus.trim_sorted_bam, irma.seg_ha_bam, irma.seg_na_bam, ""]),
            bedfile = select_first([reference_gene_locations_bed, organism_parameters.gene_locations_bed]),
            samplename = samplename,
            organism = organism_parameters.standardized_organism
        }
      }
      if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "WNV") { 
        # tasks specific to MPXV, sars-cov-2, and WNV
        call vadr_task.vadr {
          input:
            genome_fasta = select_first([consensus.consensus_seq]),
            assembly_length_unambiguous = consensus_qc.number_ATCG,
            vadr_opts = organism_parameters.vadr_opts,
            max_length = organism_parameters.vadr_maxlength
        }
      }      
      if (organism_parameters.standardized_organism == "HIV") {
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
    String read_screen_raw = raw_check_reads.read_screen
    String? read_screen_clean = clean_check_reads.read_screen
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
    String? kraken_target_organism_name = read_qc_trim.kraken_target_organism_name
    # Read QC - kraken outputs raw
    Float? kraken_human = read_qc_trim.kraken_human
    Float? kraken_sc2 = read_qc_trim.kraken_sc2
    String? kraken_target_organism = read_qc_trim.kraken_target_organism
    File? kraken_report = read_qc_trim.kraken_report
    # Read QC - kraken outputs dehosted
    Float? kraken_human_dehosted = read_qc_trim.kraken_human_dehosted
    Float? kraken_sc2_dehosted = read_qc_trim.kraken_sc2_dehosted
    String? kraken_target_organism_dehosted = read_qc_trim.kraken_target_organism_dehosted
    File? kraken_report_dehosted = read_qc_trim.kraken_report_dehosted
    # Read Alignment - Artic consensus outputs
    String assembly_fasta = select_first([consensus.consensus_seq, irma.irma_assembly_fasta, ""])
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
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([consensus.artic_pipeline_version, irma.irma_version, ""])
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
    # Nextclade outputs
    String nextclade_json = select_first([nextclade_v3.nextclade_json, ""])
    String auspice_json = select_first([ nextclade_v3.auspice_json, ""])
    String nextclade_tsv = select_first([nextclade_v3.nextclade_tsv, ""])
    String nextclade_version = select_first([nextclade_v3.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade_v3.nextclade_docker, ""])
    String nextclade_ds_tag = select_first([ha_na_nextclade_ds_tag, set_flu_ha_nextclade_values.nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag, ""])
    String nextclade_aa_subs = select_first([ha_na_nextclade_aa_subs, nextclade_output_parser.nextclade_aa_subs, ""])
    String nextclade_aa_dels = select_first([ha_na_nextclade_aa_dels, nextclade_output_parser.nextclade_aa_dels, ""])
    String nextclade_clade = select_first([nextclade_output_parser.nextclade_clade, ""])
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # Nextclade Flu outputs - NA specific columns - tamiflu mutation
    String? nextclade_tamiflu_resistance_aa_subs = nextclade_output_parser_flu_na.nextclade_tamiflu_aa_subs
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
    # Flu Outputs
    String? irma_version = irma.irma_version
    String? irma_type = irma.irma_type
    String? irma_subtype = irma.irma_subtype
    File? irma_ha_segment = irma.seg_ha_assembly
    File? irma_na_segment = irma.seg_na_assembly
    String? abricate_flu_type = abricate_flu.abricate_flu_type
    String? abricate_flu_subtype =  abricate_flu.abricate_flu_subtype
    File? abricate_flu_results = abricate_flu.abricate_flu_results
    String? abricate_flu_database =  abricate_flu.abricate_flu_database
    String? abricate_flu_version = abricate_flu.abricate_flu_version
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = flu_antiviral_substitutions.flu_A_315675_resistance
    String? flu_amantadine_resistance = flu_antiviral_substitutions.flu_amantadine_resistance
    String? flu_compound_367_resistance = flu_antiviral_substitutions.flu_compound_367_resistance
    String? flu_favipiravir_resistance = flu_antiviral_substitutions.flu_favipiravir_resistance
    String? flu_fludase_resistance = flu_antiviral_substitutions.flu_fludase_resistance
    String? flu_L_742_001_resistance = flu_antiviral_substitutions.flu_L_742_001_resistance
    String? flu_laninamivir_resistance = flu_antiviral_substitutions.flu_laninamivir_resistance
    String? flu_peramivir_resistance = flu_antiviral_substitutions.flu_peramivir_resistance
    String? flu_pimodivir_resistance = flu_antiviral_substitutions.flu_pimodivir_resistance
    String? flu_rimantadine_resistance = flu_antiviral_substitutions.flu_rimantadine_resistance
    String? flu_oseltamivir_resistance = flu_antiviral_substitutions.flu_oseltamivir_resistance
    String? flu_xofluza_resistance = flu_antiviral_substitutions.flu_xofluza_resistance
    String? flu_zanamivir_resistance = flu_antiviral_substitutions.flu_zanamivir_resistance
  }
}