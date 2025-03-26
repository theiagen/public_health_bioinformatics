version 1.0

task checkv {
  meta {
    description: "Run CheckV on viral assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3"
    Boolean? end_to_end = false
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
    # get version
    checkv -h | grep -Po "^CheckV [^:]+" | sed -e "s/CheckV //" | tee "VERSION"

    # run CheckV referencing the CheckV DB delineated by $CHECKVDB
    if [ ~{end_to_end} ]; then
      # run CheckV end-to-end mode
      checkv end_to_end \
        ~{assembly} checkv_results/ \
        -t ~{cpu}
    else
      # run CheckV completeness 
      checkv completeness \
        ~{assembly} checkv_results/ \
        -t ~{cpu}

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