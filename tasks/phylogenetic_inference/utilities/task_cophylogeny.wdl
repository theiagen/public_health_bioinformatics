version 1.0

task cophylogeny {
  input {
    File tree1
    File tree2
    
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiaphylo:0.2.0"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # grab the theiaphylo version
    phylovalidate --version | tee VERSION

    # generate a cophylogeny plot
    Rscript /theiaphylo/theiaphylo/gen_cophylo.R ~{tree1} ~{tree2}
  >>>
  runtime {
    docker: docker 
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    String theiaphylo_version = read_string("VERSION")
    File cophylogeny_with_branch_lengths = "cophylo_with_lengths.pdf"
    File cophylogeny_branch_order = "cophylo_branch_order.pdf"
  }
}
