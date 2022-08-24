version 1.0

task lissero {
  meta {
    description: "Serogroup typing prediction for Listeria monocytogenes"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/lissero:0.4.9--py_0"
    Int? cpu = 2

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
  >>>
  output {
    File lissero_results = "~{samplename}.tsv"
    String lissero_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
