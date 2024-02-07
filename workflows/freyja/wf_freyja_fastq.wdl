version 1.0

import "../../tasks/taxon_id/task_freyja_one_sample.wdl" as freyja_task
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc
import "../../tasks/alignment/task_bwa.wdl" as align
import "../../tasks/assembly/task_ivar_primer_trim.wdl" as trim_primers
import "../../tasks/task_versioning.wdl" as versioning

workflow freyja_fastq {
  input {
    File read1
    File read2
    File primer_bed
    File reference_genome
    Int trimmomatic_minlen = 25
    String samplename
  }
  call read_qc.read_QC_trim_pe as read_QC_trim {
    input:
      samplename = samplename,
      read1  = read1,
      read2  = read2,
      trim_minlen = trimmomatic_minlen,
      workflow_series = "theiacov"
  }
  call align.bwa {
    input:
      samplename = samplename,
      reference_genome=reference_genome,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean
  }
  call trim_primers.primer_trim {
    input:
      samplename = samplename,
      primer_bed = primer_bed,
      bamfile = bwa.sorted_bam
  }
  call freyja_task.freyja_one_sample as freyja {
    input:
      primer_trimmed_bam = primer_trim.trim_sorted_bam,
      samplename = samplename,
      reference_genome = reference_genome
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String freyja_fastq_wf_version = version_capture.phb_version
    String freyja_fastq_wf_analysis_date = version_capture.date
    # Read QC - fastq_scan outputs
    Int? num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    # Read QC - fastqc outputs  
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    String? fastqc_version = read_QC_trim.fastqc_version
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    # Read QC - bbduk outputs
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    # Read QC - kraken outputs
    String? kraken_version = read_QC_trim.kraken_version
    Float? kraken_human = read_QC_trim.kraken_human
    Float? kraken_sc2 = read_QC_trim.kraken_sc2
    File? kraken_report = read_QC_trim.kraken_report
    Float? kraken_human_dehosted = read_QC_trim.kraken_human_dehosted
    Float? kraken_sc2_dehosted = read_QC_trim.kraken_sc2_dehosted
    File? kraken_report_dehosted = read_QC_trim.kraken_report_dehosted
    # Read Alignment - bwa outputs
    String bwa_version = bwa.bwa_version
    String samtools_version = bwa.sam_version
    String alignment_method = "~{bwa.bwa_version}; ~{primer_trim.ivar_version}"
    File aligned_bam = primer_trim.trim_sorted_bam
    File aligned_bai = primer_trim.trim_sorted_bai
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
    String freyja_barcode_version = freyja.freyja_barcode_version
    String freyja_metadata_version = freyja.freyja_metadata_version
    File? freyja_bootstrap_lineages = freyja.freyja_bootstrap_lineages
    File? freyja_bootstrap_lineages_pdf = freyja.freyja_bootstrap_lineages_pdf
    File? freyja_bootstrap_summary = freyja.freyja_bootstrap_summary
    File? freyja_bootstrap_summary_pdf = freyja.freyja_bootstrap_summary_pdf
  }
}