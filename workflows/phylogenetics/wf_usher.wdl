version 1.0

import "../../tasks/phylogenetic_inference/task_usher.wdl" as usher_task
import "../../tasks/task_versioning.wdl" as versioning

workflow usher_workflow {
  meta {
    description: "Place samples onto global tree"
  }
  input {
    Array[File] assembly_fasta
    String tree_name
    String organism # currently available: "sars-cov-2", "mpox"
  }
  call usher_task.usher {
    input: 
      assembly_fasta = assembly_fasta,
      tree_name = tree_name,
      organism = organism
  }
  call versioning.version_capture {
  }
  output {
    File usher_uncondensed_tree = usher.usher_uncondensed_tree
    File usher_clades = usher.usher_clades
    Array[File] usher_subtrees = usher.usher_subtrees
    String usher_version  = usher.usher_version
    String usher_protobuf_version = usher.usher_protobuf_version 
    String usher_phb_version = version_capture.phb_version
    String usher_phb_analysis_date = version_capture.date
  }
}