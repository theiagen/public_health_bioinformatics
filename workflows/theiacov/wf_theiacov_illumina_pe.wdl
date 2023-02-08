version 1.0

import "../utilities/wf_read_QC_trim_pe_theiacov.wdl" as read_qc
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/assembly/task_ivar_consensus.wdl" as consensus_task
import "../../tasks/assembly/task_ivar_primer_trim.wdl" as primer_trim_task
import "../../tasks/assembly/task_ivar_variant_call.wdl" as variant_call_task
import "../../tasks/assembly/task_irma.wdl" as irma_task
import "../../tasks/quality_control/task_assembly_metrics.wdl" as assembly_metrics
import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/task_sc2_gene_coverage.wdl" as sc2_calculation
import "../../tasks/species_typing/task_quasitools.wdl" as quasitools
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
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
    String nextclade_dataset_tag = "2022-12-14T12:00:00Z"
    String? nextclade_dataset_name
    File? reference_gff
    File? reference_genome
    Int min_depth = 100
    String organism = "sars-cov-2"
    Boolean trim_primers = true
    File? adapters
    File? phix
    String nextclade_flu_h1n1_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_h3n2_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_vic_tag = "2022-06-08T12:00:00Z"
    String nextclade_flu_yam_tag = "2022-07-27T12:00:00Z"
    Int? genome_length
  }
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      samplename = samplename,
      read1_raw = read1_raw,
      read2_raw = read2_raw,
      adapters = adapters,
      phix = phix
  }
  # assembly via bwa and ivar for non-flu data
  if (organism != "flu"){
    call bwa_task.bwa {
      input:
        samplename = samplename,
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        reference_genome = reference_genome
    }
    if (trim_primers){
      call primer_trim_task.primer_trim {
        input:
          samplename = samplename,
          primer_bed = select_first([primer_bed]),
          bamfile = bwa.sorted_bam
      }
      call assembly_metrics.stats_n_coverage as stats_n_coverage_primtrim {
      input:
        samplename = samplename,
        bamfile = primer_trim.trim_sorted_bam,
      }
    }
    call variant_call_task.variant_call {
      input:
        samplename = samplename,
        bamfile = select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam]),
        reference_gff = reference_gff,
        reference_genome = reference_genome,
        variant_min_depth = min_depth
    }
    call consensus_task.consensus {
      input:
        samplename = samplename,
        bamfile = select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam]),
        reference_genome = reference_genome,
        consensus_min_depth = min_depth
    }
    call assembly_metrics.stats_n_coverage {
      input:
        samplename = samplename,
        bamfile = bwa.sorted_bam
    }
    String assembly_method_nonflu = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
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
          nextclade_flu_h1n1_tag = nextclade_flu_h1n1_tag,
          nextclade_flu_h3n2_tag = nextclade_flu_h3n2_tag,
          nextclade_flu_vic_tag = nextclade_flu_vic_tag,
          nextclade_flu_yam_tag = nextclade_flu_yam_tag,
      }
    }
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta =  select_first([consensus.consensus_seq,irma.irma_assembly_fasta]),
      reference_genome = reference_genome,
      genome_length = genome_length
  }
  # run organism-specific typing
  if (organism == "MPXV" || organism == "sars-cov-2" || organism == "flu" && select_first([abricate_flu.run_nextclade]) ) { 
    # tasks specific to either MPXV, sars-cov-2, or flu
    call nextclade.nextclade_one_sample {
      input:
        genome_fasta = select_first([consensus.consensus_seq, irma.seg4_ha_assembly]),
        dataset_name = select_first([abricate_flu.nextclade_name, nextclade_dataset_name, organism]),
        dataset_reference = select_first([abricate_flu.nextclade_ref, nextclade_dataset_reference]),
        dataset_tag = select_first([abricate_flu.nextclade_ds_tag, nextclade_dataset_tag])
    }
    call nextclade.nextclade_output_parser_one_sample {
      input:
      nextclade_tsv = nextclade_one_sample.nextclade_tsv
    }
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = select_first([consensus.consensus_seq])
  }
    call sc2_calculation.sc2_gene_coverage {
      input: 
        samplename = samplename,
        bamfile = select_first([bwa.sorted_bam]),
        min_depth = min_depth
    }
  }
  if (organism == "MPXV" || organism == "sars-cov-2" || organism == "WNV"){ 
    # tasks specific to MPXV, sars-cov-2, and WNV
    call vadr_task.vadr {
      input:
        genome_fasta = select_first([consensus.consensus_seq]),
        assembly_length_unambiguous = consensus_qc.number_ATCG
    }
  }
  if (organism == "HIV") {
    call quasitools.quasitools_illumina_pe {
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
    File read1_dehosted = read_QC_trim.read1_dehosted
    File read2_dehosted = read_QC_trim.read2_dehosted
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    Int num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String fastq_scan_version = read_QC_trim.fastq_scan_version
    Int num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String trimmomatic_version = read_QC_trim.trimmomatic_version
    String bbduk_docker = read_QC_trim.bbduk_docker
    String kraken_version = read_QC_trim.kraken_version
    Float kraken_human = read_QC_trim.kraken_human
    Float kraken_sc2 = read_QC_trim.kraken_sc2
    String? kraken_target_org = read_QC_trim.kraken_target_org
    String? kraken_target_org_name = read_QC_trim.kraken_target_org_name
    File kraken_report = read_QC_trim.kraken_report
    Float kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    String? kraken_target_org_dehosted =read_QC_trim.kraken_target_org_dehosted
    File kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment
    String? bwa_version = bwa.bwa_version
    String? samtools_version = bwa.sam_version
    File? read1_aligned = bwa.read1_aligned
    File? read2_aligned = bwa.read2_aligned
    String assembly_method = "TheiaCoV (~{version_capture.phb_version}): " + select_first([assembly_method_nonflu,irma.irma_version])
    String aligned_bam =  select_first([primer_trim.trim_sorted_bam,bwa.sorted_bam,""]) # set default values for select_first() to avoid workflow failures
    String aligned_bai = select_first([primer_trim.trim_sorted_bai,bwa.sorted_bai,""])
    Float? primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String? ivar_version_primtrim = primer_trim.ivar_version
    String? samtools_version_primtrim = primer_trim.samtools_version
    String? primer_bed_name = primer_trim.primer_bed_name
    File? ivar_tsv = variant_call.sample_variants_tsv
    File? ivar_vcf = variant_call.sample_variants_vcf
    String? ivar_variant_version = variant_call.ivar_version
    # Assembly QC
    String assembly_fasta = select_first([consensus.consensus_seq,irma.irma_assembly_fasta,""])
    String? ivar_version_consensus = consensus.ivar_version
    String? samtools_version_consensus = consensus.samtools_version
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate =  consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage =  consensus_qc.percent_reference_coverage
    Int consensus_n_variant_min_depth = min_depth
    # Alignment QC
    File? consensus_stats = stats_n_coverage.stats
    File? consensus_flagstat = stats_n_coverage.flagstat
    String meanbaseq_trim = select_first([stats_n_coverage_primtrim.meanbaseq, stats_n_coverage.meanbaseq,""])
    String meanmapq_trim = select_first([stats_n_coverage_primtrim.meanmapq, stats_n_coverage.meanmapq,""])
    String assembly_mean_coverage = select_first([stats_n_coverage_primtrim.depth, stats_n_coverage.depth,""])
    String? samtools_version_stats = stats_n_coverage.samtools_version
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
    String nextclade_json = select_first([nextclade_one_sample.nextclade_json,""])
    String auspice_json = select_first([ nextclade_one_sample.auspice_json,""])
    String nextclade_tsv = select_first([nextclade_one_sample.nextclade_tsv,""])
    String nextclade_version = select_first([nextclade_one_sample.nextclade_version,""])
    String nextclade_docker = select_first([nextclade_one_sample.nextclade_docker,""])
    String nextclade_ds_tag = select_first([abricate_flu.nextclade_ds_tag, nextclade_dataset_tag,""])
    String nextclade_aa_subs = select_first([nextclade_output_parser_one_sample.nextclade_aa_subs,""])
    String nextclade_aa_dels = select_first([nextclade_output_parser_one_sample.nextclade_aa_dels,""])
    String nextclade_clade = select_first([nextclade_output_parser_one_sample.nextclade_clade,""])
    String? nextclade_lineage = nextclade_output_parser_one_sample.nextclade_lineage
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