version 1.0

import "../../tasks/quality_control/basic_statistics/task_fastq_scan.wdl" as fastq_scan
import "../../tasks/species_typing/mycobacterium/task_clockwork.wdl" as clockwork_task
import "../../tasks/species_typing/mycobacterium/task_tbprofiler.wdl" as tbprofiler_task
import "../../tasks/species_typing/mycobacterium/task_tbp_parser.wdl" as tbp_parser_task
import "../../tasks/quality_control/read_filtering/task_trimmomatic.wdl" as trimmomatic_task
import "../../tasks/task_versioning.wdl" as versioning

workflow tbprofiler_tngs {
  meta {
    description: "Runs trimmomatic QC, tbprofiler, and tbp-parser on tNGS TB data"
  }
  input {
    File read1
    File read2
    String samplename
    
    String? tbdb_branch

    Boolean run_trimmomatic = true
    Boolean run_clockwork = false
  }
  call versioning.version_capture {
    input:
  }
  call fastq_scan.fastq_scan_pe as fastq_scan_raw {
      input:
        read1 = read1,
        read2 = read2
  }
  if (run_trimmomatic) {
    call trimmomatic_task.trimmomatic_pe {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename
    }
    call fastq_scan.fastq_scan_pe as fastq_scan_clean {
      input:
        read1 = trimmomatic_pe.read1_trimmed,
        read2 = trimmomatic_pe.read2_trimmed
    }
  }
  if (run_clockwork) {
    call clockwork_task.clockwork_decon_reads {
      input: 
        read1 = select_first([trimmomatic_pe.read1_trimmed, read1]),
        read2 = select_first([trimmomatic_pe.read2_trimmed, read2]),
        samplename = samplename
    } 
  }
  call tbprofiler_task.tbprofiler {
    input:
      read1 = select_first([clockwork_decon_reads.clockwork_cleaned_read1, trimmomatic_pe.read1_trimmed, read1]),
      read2 = select_first([clockwork_decon_reads.clockwork_cleaned_read2, trimmomatic_pe.read2_trimmed, read2]),
      samplename = samplename,
      tbdb_branch = tbdb_branch
  }
  call tbp_parser_task.tbp_parser {
    input:
      tbprofiler_json = tbprofiler.tbprofiler_output_json,
      tbprofiler_bam = tbprofiler.tbprofiler_output_bam,
      tbprofiler_bai = tbprofiler.tbprofiler_output_bai,
      samplename = samplename,
      tngs_data = true
  }
  output {
    # fastq_scan raw (per read stats)
    Int? fastq_scan_num_reads_raw1 = fastq_scan_raw.read1_seq
    Int? fastq_scan_num_reads_raw2 = fastq_scan_raw.read2_seq
    String? fastq_scan_raw_pairs = fastq_scan_raw.read_pairs
    File? fastq_scan_raw1_json = fastq_scan_raw.read1_fastq_scan_json
    File? fastq_scan_raw2_json = fastq_scan_raw.read2_fastq_scan_json
    String? fastq_scan_version = fastq_scan_raw.version
    String? fastq_scan_docker = fastq_scan_raw.fastq_scan_docker
    # fastq_scan clean (per read stats)
    Int? fastq_scan_num_reads_clean1 = fastq_scan_clean.read1_seq
    Int? fastq_scan_num_reads_clean2 = fastq_scan_clean.read2_seq
    String? fastq_scan_clean_pairs = fastq_scan_clean.read_pairs
    File? fastq_scan_clean1_json = fastq_scan_clean.read1_fastq_scan_json
    File? fastq_scan_clean2_json = fastq_scan_clean.read2_fastq_scan_json
    # trimmomatic outputs
    File? trimmomatic_read1_trimmed = trimmomatic_pe.read1_trimmed
    File? trimmomatic_read2_trimmed = trimmomatic_pe.read2_trimmed
    File? trimmomatic_stats = trimmomatic_pe.trimmomatic_stats
    String? trimmomatic_version = trimmomatic_pe.version
    String? trimmomatic_docker = trimmomatic_pe.trimmomatic_docker
    # clockwork outputs
    File? clockwork_decontaminated_read1 = clockwork_decon_reads.clockwork_cleaned_read1
    File? clockwork_decontaminated_read2 = clockwork_decon_reads.clockwork_cleaned_read2
    String? clockwork_version = clockwork_decon_reads.clockwork_version
    # tbprofiler outputs
    File tbprofiler_report_csv = tbprofiler.tbprofiler_output_csv
    File tbprofiler_report_tsv = tbprofiler.tbprofiler_output_tsv
    File tbprofiler_report_json = tbprofiler.tbprofiler_output_json
    File tbprofiler_output_alignment_bam = tbprofiler.tbprofiler_output_bam
    File tbprofiler_output_alignment_bai = tbprofiler.tbprofiler_output_bai
    String tbprofiler_version = tbprofiler.version
    String tbprofiler_main_lineage = tbprofiler.tbprofiler_main_lineage
    String tbprofiler_sub_lineage = tbprofiler.tbprofiler_sub_lineage
    String tbprofiler_dr_type = tbprofiler.tbprofiler_dr_type
    String tbprofiler_num_dr_variants = tbprofiler.tbprofiler_num_dr_variants
    String tbprofiler_num_other_variants = tbprofiler.tbprofiler_num_other_variants
    String tbprofiler_resistance_genes = tbprofiler.tbprofiler_resistance_genes
    Float tbprofiler_median_depth = tbprofiler.tbprofiler_median_depth
    Float tbprofiler_pct_reads_mapped = tbprofiler.tbprofiler_pct_reads_mapped
    # tbp_parser outputs
    File tbp_parser_looker_report_csv = tbp_parser.tbp_parser_looker_report_csv
    File tbp_parser_laboratorian_report_csv = tbp_parser.tbp_parser_laboratorian_report_csv
    File tbp_parser_lims_report_csv = tbp_parser.tbp_parser_lims_report_csv
    File tbp_parser_coverage_report = tbp_parser.tbp_parser_coverage_report
    Float tbp_parser_genome_percent_coverage = tbp_parser.tbp_parser_genome_percent_coverage
    Float tbp_parser_average_genome_depth = tbp_parser.tbp_parser_average_genome_depth
    String tbp_parser_version = tbp_parser.tbp_parser_version
    String tbp_parser_docker = tbp_parser.tbp_parser_docker
    # version capture outputs
    String tbprofiler_tngs_wf_analysis_date = version_capture.date
    String tbprofiler_tngs_wf_version = version_capture.phb_version
  }
}