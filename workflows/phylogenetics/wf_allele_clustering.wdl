version 1.0

import "../../tasks/phylogenetic_inference/task_allele_clustering.wdl" as allele_clustering_task
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
  call allele_clustering_task.allele_clustering {
    input:
      allele_jsons = allele_jsons,
      tree_algorithm = tree_building_algorithm,
      distance_metric = distance_algorithm,
      tree_name = tree_name
  }
  output {
    File concatenated_allele_jsons = allele_clustering.concatenated_jsons
    File allele_clustering_tree = allele_clustering.tree
  }
}
