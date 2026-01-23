version 1.0

import "../../tasks/alignment/task_bwa.wdl" as align
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl" as trim_primers
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja.wdl" as freyja_task
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_pe
import "../utilities/wf_read_QC_trim_se.wdl" as read_qc_se
import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc_ont
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as task_parse_mapping
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/quality_control/basic_statistics/task_qualimap.wdl" as qualimap_task

workflow freyja_fastq {
  input {
    File read1
    File? read2
    File? primer_bed
    File reference_genome
    File? reference_gff
    Int trimmomatic_min_length = 25
    String samplename
    String freyja_pathogen = "SARS-CoV-2"
    File? freyja_barcodes
    File? freyja_lineage_metadata
    Int? depth_cutoff
    Boolean ont = false
    String kraken2_target_organism = "Severe acute respiratory syndrome coronavirus 2"
    # qc check parameters
    File? qc_check_table
    # run qualimap
    Boolean run_qualimap = true
  }
  if (defined(read2)) {
    call read_qc_pe.read_QC_trim_pe as read_QC_trim_pe {
      input:
        samplename = samplename,
        read1  = read1,
        read2  = select_first([read2]),
        trim_min_length = trimmomatic_min_length,
        workflow_series = "theiacov",
        target_organism = kraken2_target_organism
    }
  }
  if (! defined(read2) && ! ont) {
    call read_qc_se.read_QC_trim_se as read_QC_trim_se {
      input:
        samplename = samplename,
        read1  = read1,
        trim_min_length = trimmomatic_min_length,
        workflow_series = "theiacov",
        target_organism = kraken2_target_organism
    }
  }
  if (ont) {
    call fasta_utilities_task.get_fasta_genome_size {
      input:
        fasta = reference_genome
    } 
    call nanoplot_task.nanoplot as nanoplot_raw {
      input:
        read1 = read1,
        samplename = samplename,
        est_genome_length = get_fasta_genome_size.fasta_length
    }
    call read_qc_ont.read_QC_trim_ont {
      input:
        samplename = samplename,
        read1 = read1,
        workflow_series = "theiacov",
        target_organism = kraken2_target_organism
    }
    call nanoplot_task.nanoplot as nanoplot_clean {
      input:
        read1 = read_QC_trim_ont.read1_clean,
        samplename = samplename,
        est_genome_length = get_fasta_genome_size.fasta_length
      }
    call minimap2_task.minimap2 {
      input:
        samplename = samplename,
        reference = reference_genome,
        query1 = read_QC_trim_ont.read1_clean,
        output_sam = true,
        mode = "map-ont",
        long_read_flags = false
    }
    call task_parse_mapping.sam_to_sorted_bam {
      input:
        samplename = samplename,
        sam = minimap2.minimap2_out
    }
  }
  if (! ont){
    call align.bwa {
      input:
        samplename = samplename,
        reference_genome = reference_genome,
        read1 = select_first([read_QC_trim_pe.read1_clean, read_QC_trim_se.read1_clean]),
        read2 = read_QC_trim_pe.read2_clean
    }
  }
  # Called when the primer_bed file is present, primers are trimmed and trimmed bam is passed to freyja
  if (defined(primer_bed)){
    call trim_primers.primer_trim {
      input:
        samplename = samplename,
        primer_bed = select_first([primer_bed]),
        bamfile = select_first([sam_to_sorted_bam.bam, bwa.sorted_bam])
    }
  }
  call freyja_task.freyja_one_sample as freyja {
    input:
      bamfile = select_first([primer_trim.trim_sorted_bam, sam_to_sorted_bam.bam, bwa.sorted_bam]),
      freyja_pathogen = freyja_pathogen,
      freyja_barcodes = freyja_barcodes,
      freyja_lineage_metadata = freyja_lineage_metadata,
      samplename = samplename,
      reference_genome = reference_genome,
      reference_gff = reference_gff,
      depth_cutoff = depth_cutoff
  }
  if (freyja_pathogen == "SARS-CoV-2") {
    File sc2_gene_bed = "gs://theiagen-public-resources-rp/reference_data/viral/sars-cov-2/sc2_gene_locations.bed"
    call gene_coverage_task.gene_coverage {
      input:
        bamfile = select_first([primer_trim.trim_sorted_bam, sam_to_sorted_bam.bam, bwa.sorted_bam]),
        bedfile = sc2_gene_bed,
        samplename = samplename,
        organism = "sars-cov-2"
    }
  }
  if (defined(qc_check_table)) {
    call qc_check.qc_check_phb as qc_check_task {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = freyja_pathogen,
        qc_check_criteria = {
          # Read counts - handle differences between ONT (Int from nanoplot) and Illumina (Int from fastq_scan)
          "num_reads_raw1": [select_first([nanoplot_raw.num_reads, read_QC_trim_pe.fastq_scan_raw1, read_QC_trim_pe.fastqc_raw1, read_QC_trim_se.fastq_scan_raw1]), "int", ">=", "false"],
          "num_reads_raw2": [if defined(read2) then select_first([read_QC_trim_pe.fastq_scan_raw2, read_QC_trim_pe.fastqc_raw2]) else read_QC_trim_pe.fastq_scan_raw2, "int", ">=", "false"],
          "num_reads_clean1": [select_first([nanoplot_clean.num_reads, read_QC_trim_pe.fastq_scan_clean1, read_QC_trim_pe.fastqc_clean1, read_QC_trim_se.fastq_scan_clean1]), "int", ">=", "false"],
          "num_reads_clean2": [if defined(read2) then select_first([read_QC_trim_pe.fastq_scan_clean2, read_QC_trim_pe.fastqc_clean2]) else read_QC_trim_pe.fastq_scan_clean2, "int", ">=", "false"],
          # Kraken metrics - available from all three read QC workflows
          "kraken_human": [select_first([read_QC_trim_pe.kraken_human, read_QC_trim_se.kraken_human, read_QC_trim_ont.kraken_human]), "float", "<=", "false"],
          "kraken_human_dehosted": [select_first([read_QC_trim_pe.kraken_human_dehosted, read_QC_trim_se.kraken_human_dehosted, read_QC_trim_ont.kraken_human_dehosted]), "float", "<=", "false"],
          # SC2-specific gene coverage - only available when freyja_pathogen == "SARS-CoV-2"
          "sc2_s_gene_mean_coverage": [gene_coverage.sc2_s_gene_depth, "float", ">=", "false"],
          "sc2_s_gene_percent_coverage": [gene_coverage.sc2_s_gene_percent_coverage, "float", ">=", "false"]
      }
    }
  }
  if (run_qualimap) {
    call qualimap_task.qualimap as qualimap {
      input:
        samplename = samplename,
        bam_file = select_first([primer_trim.trim_sorted_bam, sam_to_sorted_bam.bam, bwa.sorted_bam])
    }
  }
  call versioning.version_capture {
    input:
  }
  # capture correct alignment method
  String alignment_method_technology = if ont then "minimap2 ~{minimap2.minimap2_version}; ~{primer_trim.ivar_version}" else "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
  output {
    # Version Capture
    String freyja_fastq_wf_version = version_capture.phb_version
    String freyja_fastq_wf_analysis_date = version_capture.date
    # Read QC - fastq_scan outputs - Illumina PE and SE
    String fastq_scan_num_reads_raw1 = select_first([read_QC_trim_pe.fastq_scan_raw1, read_QC_trim_se.fastq_scan_raw1, ""])
    Int? fastq_scan_num_reads_raw2 = read_QC_trim_pe.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim_pe.fastq_scan_raw_pairs
    String fastq_scan_raw1_json = select_first([read_QC_trim_pe.fastq_scan_raw1_json, read_QC_trim_se.fastq_scan_raw1_json, ""])
    File? fastq_scan_raw2_json = read_QC_trim_pe.fastq_scan_raw2_json
    String fastq_scan_version = select_first([read_QC_trim_pe.fastq_scan_version, read_QC_trim_se.fastq_scan_version, ""])
    String fastq_scan_num_reads_clean1 = select_first([read_QC_trim_pe.fastq_scan_clean1, read_QC_trim_se.fastq_scan_clean1, ""])
    Int? fastq_scan_num_reads_clean2 = read_QC_trim_pe.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim_pe.fastq_scan_clean_pairs
    String fastq_scan_clean1_json = select_first([read_QC_trim_pe.fastq_scan_clean1_json, read_QC_trim_se.fastq_scan_clean1_json, ""])
    File? fastq_scan_clean2_json = read_QC_trim_pe.fastq_scan_clean2_json
    # Read QC - fastqc outputs - Illumina PE and SE
    String fastqc_num_reads_raw1 = select_first([read_QC_trim_pe.fastqc_raw1, read_QC_trim_se.fastqc_raw1, ""])
    Int? fastqc_num_reads_raw2 = read_QC_trim_pe.fastqc_raw2
    String fastqc_num_reads_raw_pairs = select_first([read_QC_trim_pe.fastqc_raw_pairs, ""])
    String fastqc_raw1_html = select_first([read_QC_trim_pe.fastqc_raw1_html, read_QC_trim_se.fastqc_raw1_html, ""])
    File? fastqc_raw2_html = read_QC_trim_pe.fastqc_raw2_html
    String fastqc_version = select_first([read_QC_trim_pe.fastqc_version, read_QC_trim_se.fastqc_version, ""])
    String fastqc_docker = select_first([read_QC_trim_pe.fastqc_docker, read_QC_trim_se.fastqc_docker, ""])
    String fastqc_num_reads_clean1 = select_first([read_QC_trim_pe.fastqc_clean1, read_QC_trim_se.fastqc_clean1, ""])
    Int? fastqc_num_reads_clean2 = read_QC_trim_pe.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim_pe.fastqc_clean_pairs
    String fastqc_clean1_html = select_first([read_QC_trim_pe.fastqc_clean1_html, read_QC_trim_se.fastqc_clean1_html, ""])
    File? fastqc_clean2_html = read_QC_trim_pe.fastqc_clean2_html
    # Read QC - nanoplot outputs - ONT
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
    # Read QC - trimmomatic outputs - Illumina PE and SE
    String trimmomatic_version = select_first([read_QC_trim_pe.trimmomatic_version, read_QC_trim_se.trimmomatic_version, ""])
    String trimmomatic_docker = select_first([read_QC_trim_pe.trimmomatic_docker, read_QC_trim_se.trimmomatic_docker, ""])
    # Read QC - fastp outputs  - Illumina PE and SE
    String fastp_version = select_first([read_QC_trim_pe.fastp_version, read_QC_trim_se.fastp_version, ""])
    String fastp_docker = select_first([read_QC_trim_pe.fastp_docker, read_QC_trim_se.fastp_docker, ""])
    String fastp_html_report = select_first([read_QC_trim_pe.fastp_html_report, read_QC_trim_se.fastp_html_report, ""])
    String fastp_json_report = select_first([read_QC_trim_pe.fastp_json_report, read_QC_trim_se.fastp_json_report, ""])
    # Read QC - nanoq - ONT
    String? nanoq_version = read_QC_trim_ont.nanoq_version
    # Read QC - bbduk outputs - Illumina PE and SE
    String bbduk_docker = select_first([read_QC_trim_pe.bbduk_docker, read_QC_trim_se.bbduk_docker, ""])
    # Read QC - clean reads - all
    File read1_clean = select_first([read_QC_trim_pe.read1_clean, read_QC_trim_se.read1_clean, read_QC_trim_ont.read1_clean])
    File? read2_clean = read_QC_trim_pe.read2_clean
    # Read QC - dehosting outputs - all
    File read1_dehosted = select_first([read_QC_trim_pe.read1_dehosted, read_QC_trim_se.read1_dehosted, read_QC_trim_ont.read1_dehosted])
    File? read2_dehosted = read_QC_trim_pe.read2_dehosted
    # Read QC - kraken outputs - all
    String kraken_version = select_first([read_QC_trim_pe.kraken_version, read_QC_trim_se.kraken_version, read_QC_trim_ont.kraken_version])
    Float kraken_human = select_first([read_QC_trim_pe.kraken_human, read_QC_trim_se.kraken_human, read_QC_trim_ont.kraken_human])
    String kraken_sc2 = select_first([read_QC_trim_pe.kraken_sc2, read_QC_trim_se.kraken_sc2, read_QC_trim_ont.kraken_sc2])
    String kraken_report = select_first([read_QC_trim_pe.kraken_report, read_QC_trim_se.kraken_report, read_QC_trim_ont.kraken_report])
    Float kraken_human_dehosted = select_first([read_QC_trim_pe.kraken_human_dehosted, read_QC_trim_se.kraken_human_dehosted, read_QC_trim_ont.kraken_human_dehosted])
    String kraken_sc2_dehosted = select_first([read_QC_trim_pe.kraken_sc2_dehosted, read_QC_trim_se.kraken_sc2_dehosted, read_QC_trim_ont.kraken_sc2_dehosted])
    File kraken_report_dehosted = select_first([read_QC_trim_pe.kraken_report_dehosted, read_QC_trim_se.kraken_report_dehosted, read_QC_trim_ont.kraken_report_dehosted])
    # Read Alignment - bwa outputs
    String? bwa_version = bwa.bwa_version
    String? alignment_method = alignment_method_technology
    # Read Alignment - minimap2 outputs
    String? minimap2_version = minimap2.minimap2_version
    String? minimap2_docker = minimap2.minimap2_docker
    # Read Alignment - samtools
    String samtools_version = select_first([sam_to_sorted_bam.samtools_version, bwa.sam_version])
    # Read Alignment - bam and bai files
    String aligned_bam = select_first([sam_to_sorted_bam.bam, primer_trim.trim_sorted_bam,""])
    String aligned_bai = select_first([sam_to_sorted_bam.bai, primer_trim.trim_sorted_bai,""])
    # Read Alignment - primer trimming outputs
    Float? primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String? ivar_version_primtrim = primer_trim.ivar_version
    String? samtools_version_primtrim = primer_trim.samtools_version
    String? primer_bed_name = primer_trim.primer_bed_name
    # Freyja Analysis outputs
    String freyja_version = freyja.freyja_version
    File freyja_variants = freyja.freyja_variants
    File freyja_depths = freyja.freyja_depths
    File freyja_demixed = freyja.freyja_demixed 
    Float freyja_coverage = freyja.freyja_coverage
    String freyja_barcode_file = select_first([freyja_barcodes, freyja.freyja_sc2_barcode_file, ""])
    String freyja_lineage_metadata_file = select_first([freyja_lineage_metadata, freyja.freyja_sc2_lineage_metadata_file, ""]) 
    String freyja_barcode_version = freyja.freyja_barcode_version
    String freyja_metadata_version = freyja.freyja_metadata_version
    String? freyja_bootstrap_lineages = freyja.freyja_bootstrap_lineages
    String? freyja_bootstrap_lineages_pdf = freyja.freyja_bootstrap_lineages_pdf
    String? freyja_bootstrap_summary = freyja.freyja_bootstrap_summary
    String? freyja_bootstrap_summary_pdf = freyja.freyja_bootstrap_summary_pdf
    File freyja_demixed_parsed = freyja.freyja_demixed_parsed
    String freyja_resid = freyja.freyja_resid
    String freyja_summarized = freyja.freyja_summarized
    String freyja_lineages = freyja.freyja_lineages
    String freyja_abundances = freyja.freyja_abundances
    # Gene Coverage outputs - SARS-CoV-2 only
    Float? sc2_s_gene_mean_coverage = gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = gene_coverage.sc2_s_gene_percent_coverage
    File? est_percent_gene_coverage_tsv = gene_coverage.est_percent_gene_coverage_tsv
    # QC Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # Qualimap outputs
    String? qualimap_version = qualimap.version
    String? qualimap_docker = qualimap.qualimap_docker
    File? qualimap_reports_bundle = qualimap.qualimap_reports_bundle
    File? qualimap_coverage_plots_html = qualimap.qualimap_coverage_plots_html
  }
}