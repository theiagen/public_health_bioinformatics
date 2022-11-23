version 1.0

import "../tasks/quality_control/task_trimmomatic.wdl" as trimmomatic
import "../tasks/quality_control/task_bbduk.wdl" as bbduk
import "../tasks/quality_control/task_fastq_scan.wdl" as fastq_scan
import "../tasks/taxon_id/task_midas.wdl" as midas

workflow read_QC_trim {
  meta {
    description: "Runs basic QC (fastq_scan), trimming (Trimmomatic), and adapter removal (bbduk) on illumina PE reads"
  }

  input {
    String  samplename
    File    read1_raw
    File    read2_raw
    Int?    trimmomatic_minlen = 75
    Int?    trimmomatic_quality_trim_score = 20
    Int?    trimmomatic_window_size = 10
    Int     bbduk_mem = 8
    Boolean call_midas = false
    File?    midas_db
  }
  call trimmomatic.trimmomatic_pe {
    input:
      samplename = samplename,
      read1 = read1_raw,
      read2 = read2_raw,
      trimmomatic_minlen = trimmomatic_minlen,
      trimmomatic_quality_trim_score = trimmomatic_quality_trim_score,
      trimmomatic_window_size = trimmomatic_window_size
  }
  call bbduk.bbduk_pe {
    input:
      samplename = samplename,
      read1_trimmed = trimmomatic_pe.read1_trimmed,
      read2_trimmed = trimmomatic_pe.read2_trimmed,
      mem_size_gb = bbduk_mem
  }
  call fastq_scan.fastq_scan_pe as fastq_scan_raw {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
  }
  call fastq_scan.fastq_scan_pe as fastq_scan_clean {
    input:
      read1 = bbduk_pe.read1_clean,
      read2 = bbduk_pe.read2_clean
  }
  if (call_midas) {
    call midas.midas as midas {
      input:
        samplename = samplename,
        read1 = read1_raw,
        read2 = read2_raw,
        midas_db = midas_db
    }
  }

  output {
    File	read1_clean	=	bbduk_pe.read1_clean
    File	read2_clean	=	bbduk_pe.read2_clean
    Int	fastq_scan_raw1	=	fastq_scan_raw.read1_seq
    Int	fastq_scan_raw2	=	fastq_scan_raw.read2_seq
    String	fastq_scan_raw_pairs	=	fastq_scan_raw.read_pairs
    Int	fastq_scan_clean1	=	fastq_scan_clean.read1_seq
    Int	fastq_scan_clean2	=	fastq_scan_clean.read2_seq
    String	fastq_scan_clean_pairs	=	fastq_scan_clean.read_pairs
    String	fastq_scan_version	=	fastq_scan_raw.version
    String	bbduk_docker	=	bbduk_pe.bbduk_docker
    String	trimmomatic_version	=	trimmomatic_pe.version
    String? midas_docker = midas.midas_docker
    File? midas_report = midas.midas_report
    String? midas_primary_genus = midas.midas_primary_genus
    String? midas_secondary_genus = midas.midas_secondary_genus
    String? midas_secondary_genus_coverage = midas.midas_secondary_genus_coverage
  }
}
