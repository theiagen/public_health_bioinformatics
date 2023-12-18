version 1.0

import "../../../tasks/phylogenetic_inference/task_shared_snps.wdl" as shared_snp_task
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
  call shared_snp_task.shared_snps {
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
    File snippy_concatenated_snps = shared_snps.snippy_concatenated_snps
    File snippy_shared_snp_table = shared_snps.snippy_shared_snp_table
  }
}