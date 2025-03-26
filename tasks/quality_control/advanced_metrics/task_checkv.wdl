version 1.0

task busco {
  meta {
    description: "Run CheckV on viral assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
    # get version
    checkv -h | grep -Po "^CheckV [^:]+" | sed -e "s/CheckV //" | tee "VERSION"
 
    checkv \
      ~{assembly} checkv_results/ \
      -t ~{cpu} \

  >>>
  output {
    String checkv_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}