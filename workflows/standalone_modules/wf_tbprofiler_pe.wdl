version 1.0


import "../../tasks/species_typing/task_tbprofiler.wdl" as tbprofiler_task
import "../../tasks/task_versioning.wdl" as versioning

workflow tbprofiler_wf {
  input {
      File read1
      File read2
      String samplename
      String mapper = "bwa"
      String caller = "bcftools"
      Int min_depth = 10
      Float min_af = 0.1
      Float min_af_pred = 0.1
      Int cov_frac_threshold = 1
    }
  call tbprofiler_task.tbprofiler {
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
    String tbprofiler_wf_version = version_capture.phb_version
    String tbprofiler_wf_analysis_date = version_capture.date
    File tbprofiler_report_csv = tbprofiler.tbprofiler_output_csv
    File tbprofiler_report_tsv = tbprofiler.tbprofiler_output_tsv
    File tbprofiler_output_alignment_bam = tbprofiler.tbprofiler_output_bam
    File tbprofiler_output_alignment_bai = tbprofiler.tbprofiler_output_bai
    String tbprofiler_version = tbprofiler.version
    String tbprofiler_main_lineage = tbprofiler.tbprofiler_main_lineage
    String tbprofiler_sub_lineage = tbprofiler.tbprofiler_sub_lineage
    String tbprofiler_dr_type = tbprofiler.tbprofiler_dr_type
    String tbprofiler_num_dr_variants = tbprofiler.tbprofiler_num_dr_variants
    String tbprofiler_num_other_variants = tbprofiler.tbprofiler_num_other_variants
    String tbprofiler_resistance_genes = tbprofiler.tbprofiler_resistance_genes
    }
 }
