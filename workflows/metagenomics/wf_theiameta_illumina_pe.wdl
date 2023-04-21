version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf
import "../../tasks/quality_control/task_retrieve_mapped.wdl" as retrieve_unmapped_task
import "../../tasks/assembly/task_shovill.wdl" as shovill_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiameta_ilumina_pe {
  meta {
    description: "Reference-based consensus calling for viral metagenomic sequencing data"
  }
  input {
    File read1
    File read2
    String samplename
    File reference
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
  }
  call read_qc_wf.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1_raw = read1,
        read2_raw = read2,
        workflow_series = "theiacov",
        trim_minlen = trim_minlen,
        trim_quality_trim_score = trim_quality_trim_score,
        trim_window_size = trim_window_size
    }
    call retrieve_unmapped_task.bowtie_retrieve_mapped_pe {
      input:
        read1 = select_first([read_QC_trim.read1_dehosted, read1]),
        read2 = select_first([read_QC_trim.read2_dehosted, read2]),
        samplename = samplename,
        reference = reference
    }
    call shovill_task.shovill_pe {
      input:
        read1_cleaned = bowtie_retrieve_mapped_pe.read1_mapped,
        read2_cleaned = bowtie_retrieve_mapped_pe.read2_mapped,
        samplename = samplename,
        assembler = "megahit"
    }
    call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiameta_illumina_pe_version = version_capture.phb_version
    String theiameta_illumina_pe_analysis_date = version_capture.date
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
    # Assembly - shovill outputs 
    File? assembly_fasta = shovill_pe.assembly_fasta
    File? contigs_gfa = shovill_pe.contigs_gfa
    File? contigs_fastg = shovill_pe.contigs_fastg
    File? contigs_lastgraph = shovill_pe.contigs_lastgraph
    String? shovill_pe_version = shovill_pe.shovill_version
    }
}