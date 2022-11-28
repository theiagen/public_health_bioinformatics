version 1.0


import "../tasks/task_taxon_id.wdl" as taxon
import "../tasks/task_versioning.wdl" as versioning

workflow tbprofiler_wf {
  input {
      File read1
      File read2
      String samplename
      String? mapper = "bwa"
      String? caller = "bcftools"
      Int? min_depth = 10
      Float? min_af = 0.1
      Float? min_af_pred = 0.1
      Int? cov_frac_threshold = 1
    }
  call taxon.tbprofiler_one_sample_pe {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      mapper = mapper,
      caller = caller,
      min_depth = min_depth,
      min_af = min_af,
      min_af_pred = min_af_pred,
      cov_frac_threshold = cov_frac_threshold
    }
  call versioning.version_capture{
    input:
  }
  output {
    String tb_profiler_wf_version = version_capture.phbg_version
    String tb_profiler_wf_analysis_date = version_capture.date
    File tb_profiler_report_csv = tbprofiler_one_sample_pe.tbprofiler_output_csv
    File tb_profiler_report_tsv = tbprofiler_one_sample_pe.tbprofiler_output_tsv
    File tbprofiler_output_alignment_bam = tbprofiler_one_sample_pe.tbprofiler_output_bam
    File tbprofiler_output_alignment_bai = tbprofiler_one_sample_pe.tbprofiler_output_bai
    String tb_profiler_version = tbprofiler_one_sample_pe.version
    String tb_profiler_main_lineage = tbprofiler_one_sample_pe.tb_profiler_main_lineage
    String tb_profiler_sub_lineage = tbprofiler_one_sample_pe.tb_profiler_sub_lineage
    String tb_profiler_dr_type = tbprofiler_one_sample_pe.tb_profiler_dr_type
    String tb_profiler_num_dr_variants = tbprofiler_one_sample_pe.tb_profiler_num_dr_variants
    String tb_profiler_num_other_variants = tbprofiler_one_sample_pe.tb_profiler_num_other_variants
    String tb_profiler_resistance_genes = tbprofiler_one_sample_pe.tb_profiler_resistance_genes
    }
 }
