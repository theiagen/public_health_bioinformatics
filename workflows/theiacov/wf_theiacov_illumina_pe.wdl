version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/task_screen.wdl" as screen
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/task_quasitools.wdl" as quasitools
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/gene_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String organism = "sars-cov-2" # options: "sars-cov-2", "HIV", "WNV", "MPXV", "flu"
    File read1_raw
    File read2_raw
    # sequencing values
    String seq_method = "ILLUMINA"
    File? primer_bed
    File? adapters
    File? phix
    # reference values
    File? reference_gff
    File? reference_genome
    Int? genome_length 
    # trimming parameters
    Boolean trim_primers = true
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    # assembly parameters
    Int min_depth = 100  # the minimum depth to use for consensus and variant calling
    Float consensus_min_freq = 0.6 # minimum frequency for a variant to be called as SNP in consensus genome
    Float variant_min_freq = 0.6 # minimum frequency for a variant to be reported in ivar outputs
    # nextclade inputs
    String nextclade_docker_image = "nextstrain/nextclade:2.14.0"
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-09-21T12:00:00Z"
    String? nextclade_dataset_name
    # nextclade flu inputs
    String nextclade_flu_h1n1_ha_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_h1n1_na_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_h3n2_ha_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_h3n2_na_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_vic_ha_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_vic_na_tag = "2023-04-02T12:00:00Z"
    String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    # read screen parameters
    Int min_reads = 113 # min basepairs / 300 (which is the longest available read length of an Illumina product)
    Int min_basepairs = 34000 # 20x coverage of hepatitis delta virus
    Int min_genome_size = 1700 # size of hepatitis delta virus
    Int max_genome_size = 2673870 # size of Pandoravirus salinus + 200 kb
    Int min_coverage = 10
    Int min_proportion = 40
    Boolean skip_screen = false
    # qc check parameters
    File? qc_check_table
  }
  call screen.check_reads as raw_check_reads {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_size = min_genome_size,
      max_genome_size = max_genome_size,
      min_coverage = min_coverage,
      min_proportion = min_proportion,
      skip_screen = skip_screen,
      workflow_series = "theiacov",
      organism = organism,
      expected_genome_size = genome_length
  }
  if (raw_check_reads.read_screen == "PASS") {
    call read_qc.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1_raw = read1_raw,
        read2_raw = read2_raw,
        adapters = adapters,
        phix = phix,
        workflow_series = "theiacov",
        trim_minlen = trim_minlen,
        trim_quality_trim_score = trim_quality_trim_score,
        trim_window_size = trim_window_size
    }
    call screen.check_reads as clean_check_reads {
      input:
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_size = min_genome_size,
        max_genome_size = max_genome_size,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        skip_screen = skip_screen,
        workflow_series = "theiacov",
        organism = organism,
        expected_genome_size = genome_length
    }
    if (clean_check_reads.read_screen == "PASS") {
      # assembly via bwa and ivar for non-flu data
      if (organism != "flu"){
        call consensus_call.ivar_consensus {
          input:
            samplename = samplename,
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            reference_genome = reference_genome,
            primer_bed = primer_bed,
            reference_gff = reference_gff,
            min_depth = min_depth,
            consensus_min_freq = consensus_min_freq,
            variant_min_freq = variant_min_freq,
            trim_primers = trim_primers
        }
      }
      # assembly via irma for flu organisms
      if (organism == "flu"){
        # flu-specific tasks
        call irma_task.irma {
          input:
            read1 = read_QC_trim.read1_clean,
            read2 = read_QC_trim.read2_clean,
            samplename = samplename,
            seq_method = seq_method
        }
        if (defined(irma.irma_assemblies)) {
          call abricate.abricate_flu {
            input:
              assembly = select_first([irma.irma_assembly_fasta]),
              samplename = samplename,
              nextclade_flu_h1n1_ha_tag = nextclade_flu_h1n1_ha_tag,
              nextclade_flu_h1n1_na_tag = nextclade_flu_h1n1_na_tag,
              nextclade_flu_h3n2_ha_tag = nextclade_flu_h3n2_ha_tag,
              nextclade_flu_h3n2_na_tag = nextclade_flu_h3n2_na_tag,
              nextclade_flu_vic_ha_tag = nextclade_flu_vic_ha_tag,
              nextclade_flu_vic_na_tag = nextclade_flu_vic_na_tag,
              nextclade_flu_yam_tag = nextclade_flu_yam_tag,
          }
        }
      }
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta =  select_first([ivar_consensus.assembly_fasta,irma.irma_assembly_fasta]),
          reference_genome = reference_genome,
          genome_length = genome_length
      }
      # run organism-specific typing
      if (organism == "MPXV" || organism == "sars-cov-2" || organism == "flu" && select_first([abricate_flu.run_nextclade]) ) { 
        # tasks specific to either MPXV, sars-cov-2, or flu
        call nextclade_task.nextclade {
          input:
            docker = nextclade_docker_image,
            genome_fasta = select_first([ivar_consensus.assembly_fasta, irma.seg_ha_assembly]),
            dataset_name = select_first([abricate_flu.nextclade_name_ha, nextclade_dataset_name, organism]),
            dataset_reference = select_first([abricate_flu.nextclade_ref_ha, nextclade_dataset_reference]),
            dataset_tag = select_first([abricate_flu.nextclade_ds_tag_ha, nextclade_dataset_tag])
        }
        call nextclade_task.nextclade_output_parser {
          input:
            nextclade_tsv = nextclade.nextclade_tsv,
            organism = organism
        }
      }
      if (organism == "flu" &&  select_first([abricate_flu.run_nextclade]) && defined(irma.seg_na_assembly)) { 
        # tasks specific to flu NA - run nextclade a second time
        call nextclade_task.nextclade as nextclade_flu_na {
          input:
            docker = nextclade_docker_image,
            genome_fasta = select_first([irma.seg_na_assembly]),
            dataset_name = select_first([abricate_flu.nextclade_name_na, nextclade_dataset_name, organism]),
            dataset_reference = select_first([abricate_flu.nextclade_ref_na, nextclade_dataset_reference]),
            dataset_tag = select_first([abricate_flu.nextclade_ds_tag_na, nextclade_dataset_tag])
        }
        call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_na {
          input:
            nextclade_tsv = nextclade_flu_na.nextclade_tsv,
            organism = organism,
            NA_segment = true
        }
        # concatenate tag, aa subs and aa dels for HA and NA segments
        String ha_na_nextclade_ds_tag= "~{abricate_flu.nextclade_ds_tag_ha + ',' + abricate_flu.nextclade_ds_tag_na}"
        String ha_na_nextclade_aa_subs= "~{nextclade_output_parser.nextclade_aa_subs + ',' + nextclade_output_parser_flu_na.nextclade_aa_subs}"
        String ha_na_nextclade_aa_dels= "~{nextclade_output_parser.nextclade_aa_dels + ',' + nextclade_output_parser_flu_na.nextclade_aa_dels}"
      }
      if (organism == "sars-cov-2") {
        # sars-cov-2 specific tasks
        call pangolin.pangolin4 {
          input:
            samplename = samplename,
            fasta = select_first([ivar_consensus.assembly_fasta])
      }
        call sc2_calculation.sc2_gene_coverage {
          input: 
            samplename = samplename,
            bamfile = select_first([ivar_consensus.aligned_bam]),
            min_depth = min_depth
        }
      }
      if (organism == "MPXV" || organism == "sars-cov-2" || organism == "WNV"){ 
        # tasks specific to MPXV, sars-cov-2, and WNV
        call vadr_task.vadr {
          input:
            genome_fasta = select_first([ivar_consensus.assembly_fasta]),
            assembly_length_unambiguous = consensus_qc.number_ATCG
        }
      }
      if (organism == "HIV") {
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
            expected_taxon = organism,
            num_reads_raw1 = read_QC_trim.fastq_scan_raw1,
            num_reads_raw2 = read_QC_trim.fastq_scan_raw2,
            num_reads_clean1 = read_QC_trim.fastq_scan_clean1,
            num_reads_clean2 = read_QC_trim.fastq_scan_clean2,
            kraken_human = read_QC_trim.kraken_human,
            # kraken_sc2 = read_QC_trim.kraken_sc2,
            # kraken_target_org = read_QC_trim.kraken_target_org,
            kraken_human_dehosted = read_QC_trim.kraken_human_dehosted,
            # kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted,
            # kraken_target_org_dehosted =read_QC_trim.kraken_target_org_dehosted,
            meanbaseq_trim = ivar_consensus.meanbaseq_trim,
            assembly_mean_coverage = ivar_consensus.assembly_mean_coverage,
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
    String theiacov_illumina_pe_version = version_capture.phb_version
    String theiacov_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String  seq_platform = seq_method
    # Sample Screening
    String raw_read_screen = raw_check_reads.read_screen
    String? clean_read_screen = clean_check_reads.read_screen
    # Read QC - fastq_scan outputs
    Int? num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
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
    String? kraken_target_org = read_QC_trim.kraken_target_org
    String? kraken_target_org_name = read_QC_trim.kraken_target_org_name
    File? kraken_report = read_QC_trim.kraken_report
    Float? kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float? kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String? kraken_target_org_dehosted =read_QC_trim.kraken_target_org_dehosted
    File? kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment - bwa outputs
    String? bwa_version = ivar_consensus.bwa_version
    String? samtools_version = ivar_consensus.samtools_version
    File? read1_aligned = ivar_consensus.read1_aligned
    File? read2_aligned = ivar_consensus.read2_aligned
    String aligned_bam = select_first([ivar_consensus.aligned_bam, ""])
    String aligned_bai = select_first([ivar_consensus.aligned_bai, ""])
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
    String assembly_mean_coverage = select_first([ivar_consensus.assembly_mean_coverage, ""])
    String? samtools_version_stats = ivar_consensus.samtools_version_stats
    # Read Alignment - consensus assembly summary outputs
    Int? number_N = consensus_qc.number_N
    Int? assembly_length_unambiguous = consensus_qc.number_ATCG
    Int? number_Degenerate =  consensus_qc.number_Degenerate
    Int? number_Total = consensus_qc.number_Total
    Float? percent_reference_coverage =  consensus_qc.percent_reference_coverage
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
    String nextclade_json = select_first([nextclade.nextclade_json, ""])
    String auspice_json = select_first([ nextclade.auspice_json, ""])
    String nextclade_tsv = select_first([nextclade.nextclade_tsv, ""])
    String nextclade_version = select_first([nextclade.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade.nextclade_docker, ""])
    String nextclade_ds_tag = select_first([ha_na_nextclade_ds_tag, abricate_flu.nextclade_ds_tag_ha, nextclade_dataset_tag, ""])
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