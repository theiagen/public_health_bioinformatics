version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_phylovalidate.wdl" as phylovalidate

workflow phylocompare {
  input {
    String tree1_path # consider using NEWICK as name
    String tree2_path

    Boolean? unrooted = true
    Float? rf_max_distance
  }
  call versioning.version_capture {
    input:
  }
  call phylovalidate.phylovalidate {
    input:
        tree1_path = tree1_path,
        tree2_path = tree2_path,
        rf_max_distance = rf_max_distance,
        unrooted = unrooted
  }
  output {
    String phylovalidate_version = version_capture.phb_version
    String ete3_version = phylovalidate.ete3_version
    File phylocompare_report = phylovalidate.summary_report
    Float rf_distance = phylovalidate.rf_distance
    String validation = phylovalidate.phylovalidate
  }
}
