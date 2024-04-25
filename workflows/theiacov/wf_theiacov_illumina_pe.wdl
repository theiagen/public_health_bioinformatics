version 1.0

import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate
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
import "../utilities/wf_influenza_antiviral_substitutions.wdl" as flu_antiviral
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
    Int min_depth = 100  # the minimum depth to use for consensus and variant calling
    Float consensus_min_freq = 0.6 # minimum frequency for a variant to be called as SNP in consensus genome
    Float variant_min_freq = 0.6 # minimum frequency for a variant to be reported in ivar outputs
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_options
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
      primer_bed_file = primer_bed,
      pangolin_docker_image = pangolin_docker_image,
      kraken_target_organism_input = target_organism
  }
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
      skip_screen = skip_screen,
      workflow_series = "theiacov",
      organism = organism_parameters.standardized_organism,
      expected_genome_length = organism_parameters.genome_length
  }
  if (raw_check_reads.read_screen == "PASS") {
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
        skip_screen = skip_screen,
        workflow_series = "theiacov",
        organism = organism_parameters.standardized_organism,
        expected_genome_length = organism_parameters.genome_length
    }
    if (clean_check_reads.read_screen == "PASS") {
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
            min_depth = min_depth,
            consensus_min_freq = consensus_min_freq,
            variant_min_freq = variant_min_freq,
            trim_primers = trim_primers
        }
      }
      # assembly via irma for flu organisms
      if (organism_parameters.standardized_organism == "flu") {
        # flu-specific tasks
        call irma_task.irma {
          input:
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            samplename = samplename,
            seq_method = seq_method
        }
        # can be redone later to accomodate processing of HA and NA bams together in the task, perhaps with an organism flag
        if (defined(irma.seg_ha_bam)) {
          call assembly_metrics.stats_n_coverage as ha_assembly_coverage {
            input:
              bamfile = select_first([irma.seg_ha_bam]),
              samplename = samplename
          }
        }
        if (defined(irma.seg_na_bam)) {
          call assembly_metrics.stats_n_coverage as na_assembly_coverage {
            input:
              bamfile = select_first([irma.seg_na_bam]),
              samplename = samplename
          }
        }
        String ha_na_assembly_coverage = "HA:" + select_first([ha_assembly_coverage.depth, ""]) + ", " + "NA:" + select_first([na_assembly_coverage.depth, ""])
        if (defined(irma.irma_assemblies)) {
          call abricate.abricate_flu {
            input:
              assembly = select_first([irma.irma_assembly_fasta]),
              samplename = samplename
          }
          call set_organism_defaults.organism_parameters as set_flu_na_nextclade_values {
            input:
              organism = organism_parameters.standardized_organism,
              flu_segment = "NA",
              flu_subtype = irma.irma_subtype,
              # including these to block from terra
              reference_gff_file = reference_gff,
              reference_genome = reference_genome,
              genome_length_input = genome_length,
              nextclade_dataset_tag_input = nextclade_dataset_tag,
              nextclade_dataset_name_input = nextclade_dataset_name,
              vadr_max_length = vadr_max_length,
              vadr_skip_length = vadr_skip_length,
              vadr_options = vadr_options,
              vadr_mem = vadr_memory,
              primer_bed_file = primer_bed,
              gene_locations_bed_file = reference_gene_locations_bed,
              pangolin_docker_image = pangolin_docker_image,
              kraken_target_organism_input = target_organism,
              hiv_primer_version = "N/A"
          }
          call set_organism_defaults.organism_parameters as set_flu_ha_nextclade_values {
            input:
              organism = organism_parameters.standardized_organism,
              flu_segment = "HA",
              flu_subtype = irma.irma_subtype,
              # including these to block from terra
              reference_gff_file = reference_gff,
              reference_genome = reference_genome,
              genome_length_input = genome_length,
              nextclade_dataset_tag_input = nextclade_dataset_tag,
              nextclade_dataset_name_input = nextclade_dataset_name,     
              vadr_max_length = vadr_max_length,
              vadr_skip_length = vadr_skip_length,
              vadr_options = vadr_options,
              primer_bed_file = primer_bed,
              gene_locations_bed_file = reference_gene_locations_bed,
              pangolin_docker_image = pangolin_docker_image,
              kraken_target_organism_input = target_organism,
              hiv_primer_version = "N/A"
          }
          # these are necessary because these are optional values and cannot be directly compared in before the nextclade task. checking for variable definition can be done though, which is why we create variables here
          if (set_flu_na_nextclade_values.nextclade_dataset_tag == "NA") {
            Boolean do_not_run_flu_na_nextclade = true
          }
          if (set_flu_ha_nextclade_values.nextclade_dataset_tag == "NA") {
            Boolean do_not_run_flu_ha_nextclade = true
          }
        }       
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
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta =  select_first([ivar_consensus.assembly_fasta,irma.irma_assembly_fasta]),
          reference_genome = organism_parameters.reference,
          genome_length = organism_parameters.genome_length
      }
      # run organism-specific typing
      if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || (organism_parameters.standardized_organism == "flu" && defined(irma.seg_ha_assembly) && ! defined(do_not_run_flu_ha_nextclade))) { 
        # tasks specific to either MPXV, sars-cov-2, flu, or RSV-A/RSV-B
        call nextclade_task.nextclade_v3 {
          input:
            genome_fasta = select_first([irma.seg_ha_assembly, ivar_consensus.assembly_fasta]),
            dataset_name = select_first([set_flu_ha_nextclade_values.nextclade_dataset_name, organism_parameters.nextclade_dataset_name]),
            dataset_tag = select_first([set_flu_ha_nextclade_values.nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag])
        }
        call nextclade_task.nextclade_output_parser {
          input:
            nextclade_tsv = nextclade_v3.nextclade_tsv,
            organism = organism_parameters.standardized_organism
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
            organism = organism_parameters.standardized_organism
        }
        # concatenate tag, aa subs and aa dels for HA and NA segments
        String ha_na_nextclade_ds_tag = "~{set_flu_ha_nextclade_values.nextclade_dataset_tag + ',' + set_flu_na_nextclade_values.nextclade_dataset_tag}"
        String ha_na_nextclade_aa_subs = "~{nextclade_output_parser.nextclade_aa_subs + ',' + nextclade_output_parser_flu_na.nextclade_aa_subs}"
        String ha_na_nextclade_aa_dels = "~{nextclade_output_parser.nextclade_aa_dels + ',' + nextclade_output_parser_flu_na.nextclade_aa_dels}"
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
            bamfile = select_first([ivar_consensus.aligned_bam, irma.seg_ha_bam, irma.seg_na_bam, ""]),
            bedfile = select_first([reference_gene_locations_bed, organism_parameters.gene_locations_bed]),
            samplename = samplename,
            organism = organism_parameters.standardized_organism
        }
      }
      if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b"){ 
        # tasks specific to MPXV, sars-cov-2, WNV, flu rsv_a and rsv_b
        call vadr_task.vadr {
          input:
            genome_fasta = select_first([ivar_consensus.assembly_fasta, irma.irma_assembly_fasta]),
            assembly_length_unambiguous = consensus_qc.number_ATCG,
            vadr_opts = organism_parameters.vadr_opts,
            max_length = organism_parameters.vadr_maxlength,
            skip_length = organism_parameters.vadr_skiplength,
            memory = organism_parameters.vadr_memory
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
    String  seq_platform = seq_method
    # Sample Screening
    String read_screen_raw = raw_check_reads.read_screen
    String? read_screen_clean = clean_check_reads.read_screen
    # Read QC - fastq_scan outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
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
    Float? kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_organism = read_QC_trim.kraken_target_organism
    String? kraken_target_organism_name = read_QC_trim.kraken_target_organism_name
    File? kraken_report = read_QC_trim.kraken_report
    Float? kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float? kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
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
    # Read Alignment - assembly outputs
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([ivar_consensus.assembly_method_nonflu, irma.irma_version, ""])
    String assembly_fasta = select_first([ivar_consensus.assembly_fasta, irma.irma_assembly_fasta, ""])
    String? ivar_version_consensus = ivar_consensus.ivar_version_consensus
    String? samtools_version_consensus = ivar_consensus.samtools_version_consensus
    # Read Alignment - consensus assembly qc outputs
    Int consensus_n_variant_min_depth = min_depth
    File? consensus_stats = ivar_consensus.consensus_stats
    File? consensus_flagstat = ivar_consensus.consensus_flagstat
    String meanbaseq_trim = select_first([ivar_consensus.meanbaseq_trim, ""])
    String meanmapq_trim = select_first([ivar_consensus.meanmapq_trim, ""])
    String assembly_mean_coverage = select_first([ivar_consensus.assembly_mean_coverage, ha_na_assembly_coverage , ""])
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
    # Nextclade outputs
    String nextclade_json = select_first([nextclade_v3.nextclade_json, ""])
    String nextclade_json_flu_na = select_first([nextclade_flu_na.nextclade_json, ""])
    String auspice_json = select_first([ nextclade_v3.auspice_json, ""])
    String auspice_json_flu_na = select_first([nextclade_flu_na.auspice_json, ""])
    String nextclade_tsv = select_first([nextclade_v3.nextclade_tsv, ""])
    String nextclade_tsv_flu_na = select_first([nextclade_flu_na.nextclade_tsv, ""])
    String nextclade_version = select_first([nextclade_v3.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade_v3.nextclade_docker, ""])
    String nextclade_ds_tag = select_first([ha_na_nextclade_ds_tag, set_flu_ha_nextclade_values.nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag, ""])
    String nextclade_aa_subs = select_first([ha_na_nextclade_aa_subs, nextclade_output_parser.nextclade_aa_subs, ""])
    String nextclade_aa_dels = select_first([ha_na_nextclade_aa_dels, nextclade_output_parser.nextclade_aa_dels, ""])
    String nextclade_clade = select_first([nextclade_output_parser.nextclade_clade, ""])
    String nextclade_clade_flu_na = select_first([nextclade_output_parser_flu_na.nextclade_clade, ""])
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    String? nextclade_qc_flu_na = nextclade_output_parser_flu_na.nextclade_qc
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # Flu IRMA and Abricate Outputs
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
 
  }
}
