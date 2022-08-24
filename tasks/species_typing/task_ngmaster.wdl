version 1.0

task ngmaster {
  meta {
    description: "Multi-antigen sequence typing for Neisseria gonorrhoeae"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/ngmaster:0.5.8--pyhdfd78af_1"
    Int? cpu = 2
  }
  command <<<
    echo $(ngmaster --version 2>&1) | sed 's/^.*ngmaster //' | tee VERSION
    ngmaster \
      ~{assembly} \
      > ~{samplename}.tsv
  >>>
  output {
    File ngmaster_results = "~{samplename}.tsv"
    String ngmaster_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
