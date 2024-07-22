version 1.0

import "../../tasks/alignment/task_bwa.wdl" as align
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl" as trim_primers
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/freyja/task_freyja.wdl" as freyja_task
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_pe
import "../utilities/wf_read_QC_trim_se.wdl" as read_qc_se
import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc_ont
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as task_parse_mapping

workflow freyja_fastq {
  input {
    File read1
    File? read2
    File primer_bed
    File reference_genome
    Int trimmomatic_min_length = 25
    String samplename
    Int? depth_cutoff
    Boolean ont = false
  }
  if (defined(read2)) {
    call read_qc_pe.read_QC_trim_pe as read_QC_trim_pe {
      input:
        samplename = samplename,
        read1  = read1,
        read2  = select_first([read2]),
        trim_min_length = trimmomatic_min_length,
        workflow_series = "theiacov"
    }
  }
  if (! defined(read2) && ! ont) {
    call read_qc_se.read_QC_trim_se as read_QC_trim_se {
      input:
        samplename = samplename,
        read1  = read1,
        trim_min_length = trimmomatic_min_length,
        workflow_series = "theiacov"
    }
  }
  if (ont) {
    call read_qc_ont.read_QC_trim_ont {
      input:
        samplename = samplename,
        read1 = read1
    }
    call minimap2_task.minimap2 {
      input:
        samplename = samplename,
        reference = reference_genome,
        query1 = read_QC_trim_ont.read1_clean,
        output_sam = true,
        mode = "map-ont"
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
        read2 = select_first([read_QC_trim_pe.read2_clean])
    }
  }
  call trim_primers.primer_trim {
    input:
      samplename = samplename,
      primer_bed = primer_bed,
      bamfile = select_first([sam_to_sorted_bam.bam,bwa.sorted_bam])
  }
  call freyja_task.freyja_one_sample as freyja {
    input:
      primer_trimmed_bam = primer_trim.trim_sorted_bam,
      samplename = samplename,
      reference_genome = reference_genome,
      depth_cutoff = depth_cutoff
  }
  call versioning.version_capture {
    input:
  }
  output {
    # Version Capture
    String freyja_fastq_wf_version = version_capture.phb_version
    String freyja_fastq_wf_analysis_date = version_capture.date
    # Read QC - fastq_scan outputs - Illumina PE and SE
    String fastq_scan_num_reads_raw1 = select_first([read_QC_trim_pe.fastq_scan_raw1, read_QC_trim_se.fastq_scan_raw1, ""])
    Int? fastq_scan_num_reads_raw2 = read_QC_trim_pe.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim_pe.fastq_scan_raw_pairs
    String fastq_scan_version = select_first([read_QC_trim_pe.fastq_scan_version, read_QC_trim_se.fastq_scan_version, ""])
    String fastq_scan_num_reads_clean1 = select_first([read_QC_trim_pe.fastq_scan_clean1, read_QC_trim_se.fastq_scan_clean1, ""])
    Int? fastq_scan_num_reads_clean2 = read_QC_trim_pe.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim_pe.fastq_scan_clean_pairs
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
    # Read QC - trimmomatic outputs - Illumina PE and SE
    String trimmomatic_version = select_first([read_QC_trim_pe.trimmomatic_version, read_QC_trim_se.trimmomatic_version, ""])
    String trimmomatic_docker = select_first([read_QC_trim_pe.trimmomatic_docker, read_QC_trim_se.trimmomatic_docker, ""])
    # Read QC - fastp outputs  - Illumina PE and SE
    String fastp_version = select_first([read_QC_trim_pe.fastp_version, read_QC_trim_se.fastp_version, ""])
    String fastp_html_report = select_first([read_QC_trim_pe.fastp_html_report, read_QC_trim_se.fastp_html_report, ""])
    # Read QC - nanoq - ONT
    String? nanoq_version = read_QC_trim_ont.nanoq_version
    # Read QC - bbduk outputs - Illumina PE and SE
    String bbduk_docker = select_first([read_QC_trim_pe.bbduk_docker,read_QC_trim_se.bbduk_docker])
    # Read QC - clean reads - all
    File read1_clean = select_first([read_QC_trim_pe.read1_clean, read_QC_trim_se.read1_clean, read_QC_trim_ont.read1_clean])
    File? read2_clean = read_QC_trim_pe.read2_clean
    # Read QC - dehosting outputs - all
    File read1_dehosted = select_first([read_QC_trim_pe.read1_dehosted, read_QC_trim_se.read1_dehosted, read_QC_trim_ont.read1_dehosted])
    File? read2_dehosted = read_QC_trim_pe.read2_dehosted
    # Read QC - kraken outputs - all
    String kraken_version = select_first([read_QC_trim_pe.kraken_version, read_QC_trim_se.kraken_version, read_QC_trim_ont.kraken_version])
    Float kraken_human = select_first([read_QC_trim_pe.kraken_human, read_QC_trim_se.kraken_human, read_QC_trim_ont.kraken_human])
    Float kraken_sc2 = select_first([read_QC_trim_pe.kraken_sc2, read_QC_trim_se.kraken_sc2, read_QC_trim_ont.kraken_sc2])
    String kraken_report = select_first([read_QC_trim_pe.kraken_report, read_QC_trim_se.kraken_report, read_QC_trim_ont.kraken_report])
    Float kraken_human_dehosted = select_first([read_QC_trim_pe.kraken_human_dehosted, read_QC_trim_se.kraken_human_dehosted, read_QC_trim_ont.kraken_human_dehosted])
    Float kraken_sc2_dehosted = select_first([read_QC_trim_pe.kraken_sc2_dehosted, read_QC_trim_se.kraken_sc2_dehosted, read_QC_trim_ont.kraken_sc2_dehosted])
    File kraken_report_dehosted = select_first([read_QC_trim_pe.kraken_report_dehosted, read_QC_trim_se.kraken_report_dehosted, read_QC_trim_ont.kraken_report_dehosted])
    # Read Alignment - bwa outputs
    String? bwa_version = bwa.bwa_version
    String? alignment_method = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
    # Read Alignment - minimap2 outputs
    String? minimap2_version = minimap2.minimap2_version
    String? minimap2_docker = minimap2.minimap2_docker
    # Read Alignment - samtools
    String samtools_version = select_first([sam_to_sorted_bam.samtools_version, bwa.sam_version])
    # Read Alignment - bam and bai files
    File aligned_bam = select_first([sam_to_sorted_bam.bam, primer_trim.trim_sorted_bam])
    File aligned_bai = select_first([sam_to_sorted_bam.bai, primer_trim.trim_sorted_bai])
    # Read Alignment - primer trimming outputs
    Float primer_trimmed_read_percent = primer_trim.primer_trimmed_read_percent
    String ivar_version_primtrim = primer_trim.ivar_version
    String samtools_version_primtrim = primer_trim.samtools_version
    String primer_bed_name = primer_trim.primer_bed_name
    # Freyja Analysis outputs
    String freyja_version = freyja.freyja_version
    File freyja_variants = freyja.freyja_variants
    File freyja_depths = freyja.freyja_depths
    File freyja_demixed = freyja.freyja_demixed
    Float freyja_coverage = freyja.freyja_coverage
    File freyja_barcode_file = freyja.freyja_barcode_file
    String freyja_barcode_version = freyja.freyja_barcode_version
    String freyja_metadata_version = freyja.freyja_metadata_version
    File? freyja_bootstrap_lineages = freyja.freyja_bootstrap_lineages
    File? freyja_bootstrap_lineages_pdf = freyja.freyja_bootstrap_lineages_pdf
    File? freyja_bootstrap_summary = freyja.freyja_bootstrap_summary
    File? freyja_bootstrap_summary_pdf = freyja.freyja_bootstrap_summary_pdf
  }
}