version 1.0

import "../../../tasks/phylogenetic_inference/task_concatenate_variants.wdl" as concatenate_variants_task
import "../../../tasks/task_versioning.wdl" as versioning

workflow concatenate_variants_wf {
  meta {
    description: "Concatenate results files from Snippy variants and create a shared SNP table"
  }
  input {
    String concatenated_file_name
    Array[File] snippy_variants_results
    Array[String] samplenames

  }
  call concatenate_variants_task.concatenate_variants {
    input:
      snippy_variants_results = snippy_variants_results,
      samplenames = samplenames,
      concatenated_file_name = concatenated_file_name
    }
  call versioning.version_capture{
    input:
  }
  output {
    # version capture
    String concatenate_variants_version = version_capture.phb_version
    String concatenate_variants_analysis_date = version_capture.date

    # shared snps outputs
    File snippy_concatenated_snps = concatenate_variants.snippy_concatenated_snps
    File snippy_shared_snp_table = concatenate_variants.snippy_shared_snp_table
  }
}