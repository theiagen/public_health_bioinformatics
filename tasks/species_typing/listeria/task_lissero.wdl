version 1.0

task lissero {
  meta {
    description: "Serogroup typing prediction for Listeria monocytogenes"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/lissero:0.4.9--py_0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8

    # Parameters
    # --min_id     Minimum percent identity to accept a match [Default 95.0]
    # --min_cov    Minimum coverage of the gene to accept a match [Default 95.0]
    Float min_id = 95.0
    Float min_cov = 95.0
  }
  command <<<
    echo $(lissero --version 2>&1) | sed 's/^.*LisSero //' | tee VERSION
    lissero \
      ~{'--min_id ' + min_id} \
      ~{'--min_cov ' + min_cov} \
      ~{assembly} \
      > ~{samplename}.tsv

    # pull out serotype
    tail -n+2 ~{samplename}.tsv | cut -f2 | tee SEROTYPE
  >>>
  output {
    File lissero_results = "~{samplename}.tsv"
    String lissero_version = read_string("VERSION")
    String lissero_serotype = read_string("SEROTYPE")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
