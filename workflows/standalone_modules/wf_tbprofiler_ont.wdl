version 1.0


import "../../tasks/species_typing/task_tbprofiler.wdl" as tbprofiler
import "../../tasks/task_versioning.wdl" as versioning

workflow tbprofiler_wf {
  input {
      File reads
      String samplename
      String mapper = "minimap2"
      String caller = "bcftools"
      Int min_depth = 20
      Float min_af = 0.1
      Float min_af_pred = 0.1
      Int cov_frac_threshold = 1
    }
  call tbprofiler.tbprofiler_ont {
    input:
      reads = reads,
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
    File tbprofiler_output_alignment_bam = tbprofiler_ont.tbprofiler_output_bam
    File tbprofiler_output_alignment_bai = tbprofiler_ont.tbprofiler_output_bai
    File tbprofiler_report_csv = tbprofiler_ont.tbprofiler_output_csv
    File tbprofiler_report_tsv = tbprofiler_ont.tbprofiler_output_tsv
    String tbprofiler_version = tbprofiler_ont.version
    String tbprofiler_main_lineage = tbprofiler_ont.tbprofiler_main_lineage
    String tbprofiler_sub_lineage = tbprofiler_ont.tbprofiler_sub_lineage
    String tbprofiler_dr_type = tbprofiler_ont.tbprofiler_dr_type
    String tbprofiler_num_dr_variants = tbprofiler_ont.tbprofiler_num_dr_variants
    String tbprofiler_num_other_variants = tbprofiler_ont.tbprofiler_num_other_variants
    String tbprofiler_resistance_genes = tbprofiler_ont.tbprofiler_resistance_genes
    }
 }
