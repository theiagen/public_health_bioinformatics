version 1.0

import "../../../tasks/utilities/task_file_handling.wdl" as file_handling
import "../../../tasks/phylogenetic_inference/task_shared_variants.wdl" as shared_variants_task
import "../../../tasks/task_versioning.wdl" as versioning

workflow shared_variants_wf {
  meta {
    description: "Concatenate variants from multiple samples and create a shared variants table"
  }
  input {
    Array[File] files_to_cat
    Array[String] samplenames
    String concatenated_file_name
  }
  call file_handling.cat_files {
    input:
      files_to_cat = files_to_cat,
      samplenames = samplenames,
      concatenated_file_name = concatenated_file_name,
      concatenate_variants = true
  }
  call shared_variants_task.shared_variants {
    input:
      concatenated_variants = cat_files.concatenated_files,
      concatenated_file_name = concatenated_file_name
    }
  call versioning.version_capture{
    input:
  }
  output {
    # version capture
    String shared_variants_version = version_capture.phb_version
    String shared_variants_analysis_date = version_capture.date

    # shared snps outputs
    File concatenated_variants = cat_files.concatenated_files
    File shared_variants_table = shared_variants.shared_variants_table
  }
}