version 1.0

import "../../tasks/phylogenetic_inference/task_lyveset.wdl" as lyveset
import "../../tasks/task_versioning.wdl" as versioning

workflow lyveset_workflow {
  input {
    Array[File] read1
    Array[File] read2
    Array[String] samplename
    String dataset_name
    File reference_genome
  }
  call lyveset.lyveset {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      dataset_name = dataset_name,
      reference_genome = reference_genome
  }
  call versioning.version_capture{
    input:
  }
  output {
    String lyveset_wf_version = version_capture.phb_version
    String lyveset_wf_analysis_date = version_capture.date
    String lyveset_docker_image = lyveset.lyveset_docker_image
    File? lyveset_pairwise_matrix = lyveset.lyveset_pairwise_matrix
    File? lyveset_raxml_tree = lyveset.lyveset_raxml_tree
    File? lyveset_pooled_snps_vcf = lyveset.lyveset_pooled_snps_vcf
    File? lyveset_alignment_fasta = lyveset.lyveset_alignment_fasta 
    File lyveset_log = lyveset.lyveset_log
  }
}
