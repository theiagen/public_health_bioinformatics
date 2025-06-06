version 1.0

import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/assembly/task_semibin.wdl" as semibin_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken_task
import "../../tasks/taxon_id/contamination/task_krona.wdl" as krona_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/assembly/task_spades.wdl" as spades_task
import "../../tasks/quality_control/read_filtering/task_pilon.wdl" as pilon_task
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc_wf

workflow theiameta_illumina_pe {
  meta {
    description: "Reference-based consensus calling or de novo assembly for metagenomic sequencing data"
  }
  input {
    File read1
    File read2
    String samplename
    File? reference
    File kraken2_db = "gs://theiagen-public-resources-rp/reference_data/databases/kraken2/k2_standard_08gb_20230605.tar.gz"
    Boolean output_additional_files = false
    Int spades_memory = 32
  }
  call kraken_task.kraken2_standalone as kraken2_raw {
    input:
      samplename = samplename,
      read1 = read1,
      read2 = read2,
      kraken2_db = kraken2_db,
      kraken2_args = "",
      classified_out = "classified#.fastq",
      unclassified_out = "unclassified#.fastq"
  }
  call krona_task.krona as krona_raw {
    input:
      kraken2_report = kraken2_raw.kraken2_report,
      samplename = samplename
  }
  call read_qc_wf.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = read1,
        read2 = read2,
        workflow_series = "theiameta",
        kraken_db = kraken2_db,
        call_kraken = false,
        kraken_disk_size = 100,
        kraken_memory = 8
    }
  call kraken_task.kraken2_standalone as kraken2_clean {
    input:
      samplename = samplename,
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      kraken2_db = kraken2_db,
      kraken2_args = "",
      classified_out = "classified#.fastq",
      unclassified_out = "unclassified#.fastq"
  }
  call krona_task.krona as krona_clean {
    input:
      kraken2_report = kraken2_clean.kraken2_report,
      samplename = samplename
  }
  call spades_task.spades {
    input:
      read1 = read_QC_trim.read1_clean,
      read2 = read_QC_trim.read2_clean,
      samplename = samplename,
      spades_type = 'meta',
      memory = spades_memory
  }
  call minimap2_task.minimap2 as minimap2_assembly_correction {
    input:
      query1 = read_QC_trim.read1_clean,
      query2 = read_QC_trim.read2_clean, 
      reference = select_first([spades.assembly_fasta]),
      samplename = samplename,
      mode = "sr",
      output_sam = true
  }
  call parse_mapping_task.sam_to_sorted_bam as sort_bam_assembly_correction {
    input:
      sam = minimap2_assembly_correction.minimap2_out,
      samplename = samplename
  }
  call pilon_task.pilon {
    input:
      assembly = select_first([spades.assembly_fasta]),
      bam = sort_bam_assembly_correction.bam,
      bai = sort_bam_assembly_correction.bai,
      samplename = samplename
  }
    # if reference is provided, perform mapping of assembled contigs to 
    # reference with minimap2, and extract those as final assembly
    if (defined(reference)) {
      call minimap2_task.minimap2 as minimap2_assembly {
        input:
          query1 = pilon.assembly_fasta,
          reference = select_first([reference]),
          samplename = samplename,
          mode = "asm20",
          output_sam = false,
          long_read_flags = false
      }
      call parse_mapping_task.retrieve_aligned_contig_paf {
        input:
          paf = minimap2_assembly.minimap2_out,
          assembly = pilon.assembly_fasta,
          samplename = samplename
      }
      call parse_mapping_task.calculate_coverage_paf {
        input:
          paf = minimap2_assembly.minimap2_out
      }
    }
    call quast_task.quast {
      input:
        assembly = select_first([retrieve_aligned_contig_paf.final_assembly, pilon.assembly_fasta]),
        samplename = samplename,
        min_contig_length = 1
      }
    if (output_additional_files) {
      call minimap2_task.minimap2 as minimap2_reads {
        input:
          query1 = read_QC_trim.read1_clean,
          query2 = read_QC_trim.read2_clean, 
          reference = select_first([retrieve_aligned_contig_paf.final_assembly, pilon.assembly_fasta]),
          samplename = samplename,
          mode = "sr",
          output_sam = true,
          long_read_flags = false
      }
      call parse_mapping_task.sam_to_sorted_bam {
        input:
          sam = minimap2_reads.minimap2_out,
          samplename = samplename
      }
      call parse_mapping_task.calculate_coverage {
        input:
          bam = sam_to_sorted_bam.bam,
          bai = sam_to_sorted_bam.bai
      }
      call parse_mapping_task.retrieve_pe_reads_bam as retrieve_unaligned_pe_reads_sam {
        input:
          bam = sam_to_sorted_bam.bam,
          samplename = samplename,
          prefix = "unassembled",
          sam_flag = 4
      }
      call parse_mapping_task.retrieve_pe_reads_bam as retrieve_aligned_pe_reads_sam {
        input:
          bam = sam_to_sorted_bam.bam,
          samplename = samplename,
          sam_flag = 2,
          prefix = "assembled"
      }
      call parse_mapping_task.assembled_reads_percent {
        input:
          bam = sam_to_sorted_bam.bam,
      } 
    }
    if (! defined(reference)) {
      call bwa_task.bwa as bwa {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          reference_genome = pilon.assembly_fasta,
          samplename = samplename
      }
      call semibin_task.semibin as semibin {
        input:
          sorted_bam = bwa.sorted_bam,
          sorted_bai = bwa.sorted_bai,
          assembly_fasta = pilon.assembly_fasta,
          samplename = samplename
      }
    }
    call versioning.version_capture {
      input:
  }
  output {
    # Version capture
    String theiameta_illumina_pe_version = version_capture.phb_version
    String theiameta_illumina_pe_analysis_date = version_capture.date
    # Kraken2 outputs
    String kraken2_version = kraken2_raw.kraken2_version
    String kraken2_docker = kraken2_raw.kraken2_docker
    File kraken2_report_raw = kraken2_raw.kraken2_report
    Float kraken2_percent_human_raw = kraken2_raw.kraken2_percent_human
    File kraken2_report_clean = kraken2_clean.kraken2_report
    Float kraken2_percent_human_clean = kraken2_clean.kraken2_percent_human
    # Krona outputs
    String krona_version = krona_raw.krona_version
    String krona_docker = krona_raw.krona_docker
    File krona_html_raw = krona_raw.krona_html
    File krona_html_clean = krona_clean.krona_html
    # Read QC - dehosting outputs
    File? read1_dehosted = read_QC_trim.read1_dehosted
    File? read2_dehosted = read_QC_trim.read2_dehosted
    String? ncbi_scrub_docker = read_QC_trim.ncbi_scrub_docker
    # Read QC - fastq_scan outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? fastq_scan_docker = read_QC_trim.fastq_scan_docker
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # Read QC - bbduk outputs
    File read1_clean = read_QC_trim.read1_clean
    File read2_clean = read_QC_trim.read2_clean
    String bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - Read stats
    Float? average_read_length = read_QC_trim.average_read_length
    # MIDAS outputs
    String? midas_primary_genus = read_QC_trim.midas_primary_genus
    File? midas_report = read_QC_trim.midas_report
    # Assembly - metaspades 
    File assembly_fasta = select_first([retrieve_aligned_contig_paf.final_assembly, pilon.assembly_fasta])
    String metaspades_version = spades.spades_version
    String metaspades_docker = spades.spades_docker
    # Assembly - minimap2
    String minimap2_version = minimap2_assembly_correction.minimap2_version
    String minimap2_docker = minimap2_assembly_correction.minimap2_docker
    # Assembly - samtools
    String samtools_version = sort_bam_assembly_correction.samtools_version
    String samtools_docker = sort_bam_assembly_correction.samtools_docker
    # Assembly - pilon
    String pilon_version = pilon.pilon_version
    String pilon_docker = pilon.pilon_docker
    # Assembly QC - quast
    Int assembly_length = quast.genome_length
    Int contig_number = quast.number_contigs
    Int largest_contig = quast.largest_contig
    String quast_version = quast.version
    String quast_docker = quast.quast_docker
    # Assembly QC - minimap2
    Float? percent_coverage = calculate_coverage_paf.percent_coverage
    # Assembly QC - bedtools
    Float? assembly_mean_coverage = calculate_coverage.mean_depth_coverage
    String? bedtools_version = calculate_coverage.bedtools_version
    String? bedtools_docker = calculate_coverage.bedtools_docker
    # Read retrieval
    File? read1_unmapped = retrieve_unaligned_pe_reads_sam.read1
    File? read2_unmapped = retrieve_unaligned_pe_reads_sam.read2
    File? read1_mapped = retrieve_aligned_pe_reads_sam.read1 
    File? read2_mapped = retrieve_aligned_pe_reads_sam.read2
    # Assembly stats
    Float? percentage_mapped_reads = assembled_reads_percent.percentage_mapped
    # Binning
    String? semibin_version = semibin.semibin_version
    String? semibin_docker = semibin.semibin_docker
    Array[File]? semibin_bins = semibin.semibin_bins
    }
}
