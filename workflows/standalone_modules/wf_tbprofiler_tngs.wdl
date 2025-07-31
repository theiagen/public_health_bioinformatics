version 1.0

import "../../tasks/quality_control/read_filtering/task_trimmomatic.wdl" as trimmomatic_task
import "../../tasks/species_typing/mycobacterium/task_tbprofiler.wdl" as tbprofiler_task
import "../../tasks/species_typing/mycobacterium/task_tbp_parser.wdl" as tbp_parser_task
import "../../tasks/task_versioning.wdl" as versioning

workflow tbprofiler_tngs {
  meta {
    description: "Runs trimmomatic QC, tbprofiler, and tbp-parser on tNGS TB data"
  }
  input {
    File read1
    File? read2
    String samplename
    Int bases_to_crop = 0
    Boolean skip_trimmomatic = false
    Boolean ont_data = false
  }
  call versioning.version_capture {
    input:
  }
  if (! skip_trimmomatic && ! ont_data) {
    call trimmomatic_task.trimmomatic_pe {
      input:
        read1 = read1,
        read2 = select_first([read2]),
        samplename = samplename,
        trimmomatic_base_crop = bases_to_crop
    }
  }
  if (ont_data) {
    call tbprofiler_task.tbprofiler as tbprofiler_ont {
      input:
        read1 = read1,
        samplename = samplename
    }
  }
  if (! ont_data) {
    call tbprofiler_task.tbprofiler as tbprofiler_illumina {
      input:
        read1 = select_first([trimmomatic_pe.read1_trimmed, read1]),
        read2 = select_first([trimmomatic_pe.read2_trimmed, read2]),
        samplename = samplename
    }
  }
  call tbp_parser_task.tbp_parser {
    input:
      tbprofiler_json = select_first([tbprofiler_illumina.tbprofiler_output_json, tbprofiler_ont.tbprofiler_output_json]),
      tbprofiler_bam = select_first([tbprofiler_illumina.tbprofiler_output_bam, tbprofiler_ont.tbprofiler_output_bam]),
      tbprofiler_bai = select_first([tbprofiler_illumina.tbprofiler_output_bai, tbprofiler_ont.tbprofiler_output_bai]),
      samplename = samplename,
      tngs_data = true
  }
  output {
    # trimmomatic outputs
    File? trimmomatic_read1_trimmed = trimmomatic_pe.read1_trimmed
    File? trimmomatic_read2_trimmed = trimmomatic_pe.read2_trimmed
    File? trimmomatic_stats = trimmomatic_pe.trimmomatic_stats
    String? trimmomatic_version = trimmomatic_pe.version
    String? trimmomatic_docker = trimmomatic_pe.trimmomatic_docker
    # tbprofiler outputs
    File tbprofiler_report_csv = select_first([tbprofiler_illumina.tbprofiler_output_csv, tbprofiler_ont.tbprofiler_output_csv])
    File tbprofiler_report_tsv = select_first([tbprofiler_illumina.tbprofiler_output_tsv, tbprofiler_ont.tbprofiler_output_tsv])
    File tbprofiler_report_json = select_first([tbprofiler_illumina.tbprofiler_output_json, tbprofiler_ont.tbprofiler_output_json])
    File tbprofiler_output_alignment_bam = select_first([tbprofiler_illumina.tbprofiler_output_bam, tbprofiler_ont.tbprofiler_output_bam])
    File tbprofiler_output_alignment_bai = select_first([tbprofiler_illumina.tbprofiler_output_bai, tbprofiler_ont.tbprofiler_output_bai])
    String tbprofiler_version = select_first([tbprofiler_illumina.version, tbprofiler_ont.version])
    String tbprofiler_main_lineage = select_first([tbprofiler_illumina.tbprofiler_main_lineage, tbprofiler_ont.tbprofiler_main_lineage])
    String tbprofiler_sub_lineage = select_first([tbprofiler_illumina.tbprofiler_sub_lineage, tbprofiler_ont.tbprofiler_sub_lineage])
    String tbprofiler_dr_type = select_first([tbprofiler_illumina.tbprofiler_dr_type, tbprofiler_ont.tbprofiler_dr_type])
    String tbprofiler_num_dr_variants = select_first([tbprofiler_illumina.tbprofiler_num_dr_variants, tbprofiler_ont.tbprofiler_num_dr_variants])
    String tbprofiler_num_other_variants = select_first([tbprofiler_illumina.tbprofiler_num_other_variants, tbprofiler_ont.tbprofiler_num_other_variants])
    String tbprofiler_resistance_genes = select_first([tbprofiler_illumina.tbprofiler_resistance_genes, tbprofiler_ont.tbprofiler_resistance_genes])
    Float tbprofiler_median_depth = select_first([tbprofiler_illumina.tbprofiler_median_depth, tbprofiler_ont.tbprofiler_median_depth])
    Float tbprofiler_pct_reads_mapped = select_first([tbprofiler_illumina.tbprofiler_pct_reads_mapped, tbprofiler_ont.tbprofiler_pct_reads_mapped])
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