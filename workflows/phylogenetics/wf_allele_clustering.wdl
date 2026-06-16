version 1.0

import "../../tasks/phylogenetic_inference/task_allele_clustering.wdl" as allele_cluster
import "../../tasks/task_versioning.wdl" as versioning

workflow allele_clustering {
  meta {
    description: "PulseNet 2.0 Allele Clustering using the inputs of the PulseNet 2.0 Allele Calling algorithm"
  }
  input {
    Array[File] allele_jsons
    String tree_building_algorithm
    String distance_algorithm
    String tree_name
  }
  call versioning.version_capture {
    input:
  }
  call allele_cluster.allele_clustering as allele_clustering_task {
    input:
      allele_jsons = allele_jsons,
      tree_building_algorithm = tree_building_algorithm,
      distance_algorithm = distance_algorithm,
      tree_name = tree_name
  }
  output {
    File concatenated_allele_jsons = allele_clustering_task.concatenated_jsons
    File allele_clustering_tree = allele_clustering_task.tree
    String allele_clustering_wf_version = version_capture.phb_version
    String allele_clustering_wf_analysis_date = version_capture.date
  }
}
