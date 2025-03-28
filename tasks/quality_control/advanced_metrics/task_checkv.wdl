version 1.0

task checkv {
  meta {
    description: "Run CheckV on viral assemblies"
  }
  input {
    File assembly
    String samplename
    File checkv_db = "gs://theiagen-large-public-files-rp/terra/databases/checkv/checkv-db-v1.5.tar.gz"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
    # get version
    checkv -h | grep -Po "^CheckV [^:]+" | sed -e "s/CheckV //" | tee "VERSION"

    # extract CheckV DB
    tar -xzf ~{checkv_db}
    untarred_checkv_db=$(basename ~{checkv_db} .tar.gz)

    # run CheckV referencing the CheckV DB delineated by $CHECKVDB
    checkv end_to_end \
      ~{assembly} checkv_results/ \
      -d ${untarred_checkv_db} \
      -t ~{cpu} 

  >>>
  output {
    String checkv_version = read_string("VERSION")
    File checkv_summary = "checkv_results/quality_summary.tsv"
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