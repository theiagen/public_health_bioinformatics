version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf
import "../utilities/wf_ivar_consensus.wdl" as consensus_call
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/assembly/task_shovill.wdl" as shovill_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/task_parse_paf.wdl" as parse_paf_task
import "../../tasks/utilities/task_compare_assemblies.wdl" as compare_assemblies_task
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiameta_illumina_pe {
  meta {
    description: "Reference-based consensus calling or de novo assembly for viral metagenomic sequencing data"
  }
  input {
    File read1
    File read2
    String samplename
    File? reference
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    Int min_depth = 10  # the minimum depth to use for consensus and variant calling
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
    # if reference is provided, perform consensus assembly with ivar
    if (defined(reference)){
      call consensus_call.ivar_consensus {
        input:
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          reference_genome = reference,
          min_depth = min_depth,
          trim_primers = false
        }
      call shovill_task.shovill_pe as shovil_consensus {
        input:
          read1_cleaned = read_QC_trim.read1_clean,
          read2_cleaned = read_QC_trim.read2_clean,
          samplename = samplename,
          assembler = "megahit"
      }
      call minimap2_task.minimap2 {
        input:
          query = shovil_consensus.assembly_fasta,
          reference = reference,
          samplename = samplename
      }
      call parse_paf_task.parse_paf {
        input:
          paf = minimap2.minimap2_paf,
          assembly = shovil_consensus.assembly_fasta,
          samplename = samplename
      }
      call compare_assemblies_task.compare_assemblies {
        input:
          assembly_denovo = parse_paf.parse_paf_contigs,
          assembly_consensus = ivar_consensus.assembly_fasta,
          samplename = samplename
      }
      call consensus_qc_task.consensus_qc {
        input:
          assembly_fasta =  compare_assemblies.final_assembly,
          reference_genome = reference
      }
    }
    # otherwise, perform de novo assembly with megahit
    if (!defined(reference)) {
      call shovill_task.shovill_pe as shovil_denovo {
        input:
          read1_cleaned = read_QC_trim.read1_clean,
          read2_cleaned = read_QC_trim.read2_clean,
          samplename = samplename,
          assembler = "megahit"
      }
      call quast_task.quast {
        input:
          assembly = shovil_denovo.assembly_fasta,
          samplename = samplename
      }
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
    # Assembly - shovill/ivar outputs 
    File? assembly_fasta = select_first([compare_assemblies.final_assembly, shovil_denovo.assembly_fasta])
    String? assembly_length = select_first([consensus_qc.number_Total, quast.genome_length])
    String? shovill_pe_version = shovil_denovo.shovill_version
    Int? largest_contig = quast.largest_contig
    String? ivar_version_consensus = ivar_consensus.ivar_version_consensus
    String? samtools_version_consensus = ivar_consensus.samtools_version_consensus
    }
}