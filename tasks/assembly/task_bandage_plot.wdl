version 1.0

task bandage_plot {
    input {
      File assembly_graph_gfa
      String samplename
      Int cpu = 1
      Int memory = 4
      Int disk_size = 10
      String docker = "us-docker.pkg.dev/general-theiagen/staphb/bandage:0.8.1"
    }
    command <<< 
      set -euo pipefail
      Bandage --version | tee VERSION
      Bandage image ~{assembly_graph_gfa} ~{samplename}_bandage_plot.png
    >>>
    output {
      File plot = "~{samplename}_bandage_plot.png"
      String bandage_version = read_string("VERSION")
    }
    runtime {
      docker: "~{docker}"
      cpu: cpu
      memory: "~{memory} GB"
      disks: "local-disk " + disk_size + " HDD"
      disk: disk_size + " GB"
      maxRetries: 1
      preemptible: 0
  }
}
