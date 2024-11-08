version 1.0

import "../../tasks/utilities/file_handling/task_cat_files.wdl" as file_handling
import "../../tasks/phylogenetic_inference/utilities/task_shared_variants.wdl" as shared_variants_task
import "../../tasks/task_versioning.wdl" as versioning

workflow shared_variants_wf {
  meta {
    description: "Concatenate variants from multiple samples and create a shared variants table"
  }
  input {
    Array[File] variants_to_cat
    Array[String] samplenames
    String concatenated_file_name
  }
  String concatenated_file_name_updated = sub(concatenated_file_name, " ", "_")
  call file_handling.cat_variants {
    input:
      variants_to_cat = variants_to_cat,
      samplenames = samplenames,
      concatenated_file_name = concatenated_file_name_updated
  }
  call shared_variants_task.shared_variants {
    input:
      concatenated_variants = cat_variants.concatenated_variants,
      concatenated_file_name = concatenated_file_name_updated
    }
  call versioning.version_capture{
    input:
  }
  output {
    # version capture
    String shared_variants_version = version_capture.phb_version
    String shared_variants_analysis_date = version_capture.date

    # shared snps outputs
    File concatenated_variants = cat_variants.concatenated_variants
    File shared_variants_table = shared_variants.shared_variants_table
  }
}