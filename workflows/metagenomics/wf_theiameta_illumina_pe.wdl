version 1.0

import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf
import "../../tasks/assembly/task_megahit.wdl" as megahit_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/utilities/task_parse_paf.wdl" as parse_paf_task
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/task_versioning.wdl" as versioning

workflow theiameta_illumina_pe {
  meta {
    description: "Reference-based consensus calling or de novo assembly for metagenomic sequencing data"
  }
  input {
    File read1
    File read2
    String samplename
    File? reference
    Int trim_minlen = 75
    Int trim_quality_trim_score = 30
    Int trim_window_size = 4
    File kraken2_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/k2_standard_8gb_20210517.tar.gz"
  }
  call read_qc_wf.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1_raw = read1,
        read2_raw = read2,
        workflow_series = "theiameta",
        trim_minlen = trim_minlen,
        trim_quality_trim_score = trim_quality_trim_score,
        trim_window_size = trim_window_size,
        kraken2_db = kraken2_db
    }
    call megahit_task.megahit_pe as megahit {
      input:
        read1_cleaned = read_QC_trim.read1_clean,
        read2_cleaned = read_QC_trim.read2_clean,
        samplename = samplename
      }
    # if reference is provided, perform mapping of assembled contigs to 
    # reference with minimap2, and extract those as final assembly
    if (defined(reference)){
      call minimap2_task.minimap2 as minimap2_assembly {
        input:
          query1 = megahit.assembly_fasta,
          reference = select_first([reference]),
          samplename = samplename
      }
      call parse_paf_task.retrieve_aligned_contig_paf {
        input:
          paf = minimap2_assembly.minimap2_out,
          assembly = megahit.assembly_fasta,
          samplename = samplename
      }
      call parse_paf_task.calculate_coverage_paf {
        input:
          paf = minimap2_assembly.minimap2_out
      }
    }
    call quast_task.quast {
      input:
        assembly = select_first([retrieve_aligned_contig_paf.final_assembly, megahit.assembly_fasta]),
        samplename = samplename,
        min_contig_len = 1
      }
    call minimap2_task.minimap2 as minimap2_reads {
      input:
        query1 = read_QC_trim.read1_clean,
        query2 = read_QC_trim.read2_clean, 
        reference = select_first([retrieve_aligned_contig_paf.final_assembly, megahit.assembly_fasta]),
        samplename = samplename,
        mode = "sr",
        output_sam = true
    }
    call parse_paf_task.sam_to_sorted_bam {
      input:
        sam = minimap2_reads.minimap2_out,
        samplename = samplename
    }
    call parse_paf_task.calculate_coverage {
      input:
        bam = sam_to_sorted_bam.bam,
        bai = sam_to_sorted_bam.bai
    }
    call parse_paf_task.retrieve_pe_reads_bam as retrieve_unaligned_pe_reads_sam {
      input:
        bam = sam_to_sorted_bam.bam,
        samplename = samplename
    }
    call parse_paf_task.retrieve_pe_reads_bam as retrieve_aligned_pe_reads_sam {
      input:
        bam = sam_to_sorted_bam.bam,
        samplename = samplename,
        sam_flag = 2
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
    File? kraken_report = read_QC_trim.kraken_report
    Float? kraken_percent_human = read_QC_trim.kraken_human
    # Read QC - Read stats
    Float? average_read_length = read_QC_trim.average_read_length
    # Assembly - megahit outputs 
    File? assembly_fasta = select_first([retrieve_aligned_contig_paf.final_assembly, megahit.assembly_fasta])
    String? megahit_pe_version = megahit.megahit_version
    # Assembly QC
    Int? assembly_length = quast.genome_length
    Int? contig_number = quast.number_contigs
    Int? largest_contig = quast.largest_contig
    Float? percent_coverage = calculate_coverage_paf.percent_coverage
    Float? assembly_mean_coverage = calculate_coverage.mean_depth_coverage
    # Read retrieval
    File? read1_unmapped = retrieve_unaligned_pe_reads_sam.read1
    File? read2_unmapped = retrieve_unaligned_pe_reads_sam.read2
    File? read1_mapped = retrieve_aligned_pe_reads_sam.read1 
    File? read2_mapped = retrieve_aligned_pe_reads_sam.read2
    }
}