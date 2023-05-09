version 1.0

task minimap2 {
  meta {
    description: "Align a query genome to a reference with minimap2"
  }
  input {
    File query
    File reference
    String samplename
    String docker = "staphb/minimap2:2.22" # newer versions seem to be bugged (infinite loop)
    String mode = "asm5"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    minimap2 --version | tee VERSION

    # run minimap2
    minimap2 \
      -ax "~{mode}" \
      -t "~{cpu}" \
      "~{reference}" \
      "~{query}" > "~{samplename}"_minimap2.paf

  >>>
  output {
    File minimap2_paf = "~{samplename}_minimap2.paf"
    String minimap2_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}