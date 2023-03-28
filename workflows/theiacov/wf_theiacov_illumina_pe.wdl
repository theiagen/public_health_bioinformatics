version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/task_quasitools.wdl" as quasitools
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/gene_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_illumina_pe {
  meta {
    description: "Reference-based consensus calling for viral amplicon sequencing data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    File? primer_bed
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-02-25T12:00:00Z"
    String? nextclade_dataset_name
    File? reference_gff
    File? reference_genome
    Int min_depth = 100
    String organism = "sars-cov-2"
    Boolean trim_primers = true
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    File? adapters
    File? phix
    String nextclade_flu_h1n1_ha_tag = "2023-01-27T12:00:00Z"
    String nextclade_flu_h1n1_na_tag = "2023-01-27T12:00:00Z"
    String nextclade_flu_h3n2_ha_tag = "2023-02-01T12:00:00Z"
    String nextclade_flu_h3n2_na_tag = "2023-01-27T12:00:00Z"
    String nextclade_flu_vic_ha_tag = "2023-02-01T12:00:00Z"
    String nextclade_flu_vic_na_tag = "2023-01-27T12:00:00Z"
    String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    Int? genome_length
  }
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1_raw,
      read2_raw = read2_raw,
      trim_minlen = trim_minlen,
      trim_quality_trim_score = trim_quality_trim_score,
      trim_window_size = trim_window_size,
      adapters = adapters,
      phix = phix,
      workflow_series = "theiacov"
  }
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
    call nextclade.nextclade_one_sample {
      input:
        genome_fasta = select_first([ivar_consensus.assembly_fasta, irma.seg_ha_assembly]),
        dataset_name = select_first([abricate_flu.nextclade_name_ha, nextclade_dataset_name, organism]),
        dataset_reference = select_first([abricate_flu.nextclade_ref_ha, nextclade_dataset_reference]),
        dataset_tag = select_first([abricate_flu.nextclade_ds_tag_ha, nextclade_dataset_tag])
    }
    call nextclade.nextclade_output_parser_one_sample {
      input:
        nextclade_tsv = nextclade_one_sample.nextclade_tsv,
        organism = organism
    }
  }
  if (organism == "flu" &&  select_first([abricate_flu.run_nextclade]) && defined(irma.seg_na_assembly)) { 
    # tasks specific to flu NA - run nextclade a second time
    call nextclade.nextclade_one_sample as nextclade_one_sample_flu_na {
      input:
        genome_fasta = select_first([irma.seg_na_assembly]),
        dataset_name = select_first([abricate_flu.nextclade_name_na, nextclade_dataset_name, organism]),
        dataset_reference = select_first([abricate_flu.nextclade_ref_na, nextclade_dataset_reference]),
        dataset_tag = select_first([abricate_flu.nextclade_ds_tag_na, nextclade_dataset_tag])
    }
    call nextclade.nextclade_output_parser_one_sample as nextclade_output_parser_one_sample_flu_na {
      input:
        nextclade_tsv = nextclade_one_sample_flu_na.nextclade_tsv,
        organism = organism,
        NA_segment = true
    }
    # concatenate tag, aa subs and aa dels for HA and NA segments
    String ha_na_nextclade_ds_tag= "~{abricate_flu.nextclade_ds_tag_ha + ',' + abricate_flu.nextclade_ds_tag_na}"
    String ha_na_nextclade_aa_subs= "~{nextclade_output_parser_one_sample.nextclade_aa_subs + ',' + nextclade_output_parser_one_sample_flu_na.nextclade_aa_subs}"
    String ha_na_nextclade_aa_dels= "~{nextclade_output_parser_one_sample.nextclade_aa_dels + ',' + nextclade_output_parser_one_sample_flu_na.nextclade_aa_dels}"
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
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_illumina_pe_version = version_capture.phb_version
    String theiacov_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String  seq_platform = seq_method
    # Read QC
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker
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
    # Read Alignment
    String? bwa_version = ivar_consensus.bwa_version
    String? samtools_version = ivar_consensus.samtools_version
    File? read1_aligned = ivar_consensus.read1_aligned
    File? read2_aligned = ivar_consensus.read2_aligned
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([ivar_consensus.assembly_method_nonflu, irma.irma_version])
    String aligned_bam = select_first([ivar_consensus.aligned_bam, ""])
    String aligned_bai = select_first([ivar_consensus.aligned_bai, ""])
    Float? primer_trimmed_read_percent = ivar_consensus.primer_trimmed_read_percent
    String? ivar_version_primtrim = ivar_consensus.ivar_version_primtrim
    String? samtools_version_primtrim = ivar_consensus.samtools_version
    String? primer_bed_name = ivar_consensus.primer_bed_name
    File? ivar_tsv = ivar_consensus.ivar_tsv
    File? ivar_vcf = ivar_consensus.ivar_vcf
    String? ivar_variant_version = ivar_consensus.ivar_variant_version
    # Assembly QC
    String assembly_fasta = select_first([ivar_consensus.assembly_fasta, irma.irma_assembly_fasta, ""])
    String? ivar_version_consensus = ivar_consensus.ivar_version_consensus
    String? samtools_version_consensus = ivar_consensus.samtools_version_consensus
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate =  consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage =  consensus_qc.percent_reference_coverage
    Int consensus_n_variant_min_depth = min_depth
    # Alignment QC
    File? consensus_stats = ivar_consensus.consensus_stats
    File? consensus_flagstat = ivar_consensus.consensus_flagstat
    String meanbaseq_trim = select_first([ivar_consensus.meanbaseq_trim, ""])
    String meanmapq_trim = select_first([ivar_consensus.meanmapq_trim, ""])
    String assembly_mean_coverage = select_first([ivar_consensus.assembly_mean_coverage, ""])
    String? samtools_version_stats = ivar_consensus.samtools_version_stats
    # SC2 specific
    Float? sc2_s_gene_mean_coverage = sc2_gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = sc2_gene_coverage.sc2_s_gene_percent_coverage
    File? sc2_all_genes_percent_coverage = sc2_gene_coverage.sc2_all_genes_percent_coverage
    # SC2 Lineage Assignment
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Clade Assigment
    String nextclade_json = select_first([nextclade_one_sample.nextclade_json, ""])
    String auspice_json = select_first([ nextclade_one_sample.auspice_json, ""])
    String nextclade_tsv = select_first([nextclade_one_sample.nextclade_tsv, ""])
    String nextclade_version = select_first([nextclade_one_sample.nextclade_version, ""])
    String nextclade_docker = select_first([nextclade_one_sample.nextclade_docker, ""])
    String nextclade_ds_tag = select_first([ha_na_nextclade_ds_tag, abricate_flu.nextclade_ds_tag_ha, nextclade_dataset_tag, ""])
    String nextclade_aa_subs = select_first([ha_na_nextclade_aa_subs, nextclade_output_parser_one_sample.nextclade_aa_subs, ""])
    String nextclade_aa_dels = select_first([ha_na_nextclade_aa_dels, nextclade_output_parser_one_sample.nextclade_aa_dels, ""])
    String nextclade_clade = select_first([nextclade_output_parser_one_sample.nextclade_clade, ""])
    String? nextclade_lineage = nextclade_output_parser_one_sample.nextclade_lineage
    # NA specific columns - tamiflu mutation
    String? nextclade_tamiflu_resistance_aa_subs = nextclade_output_parser_one_sample_flu_na.nextclade_tamiflu_aa_subs
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # Flu Outputs
    String? irma_version = irma.irma_version
    String? irma_type = irma.irma_type
    String? irma_subtype = irma.irma_subtype
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
  }
}