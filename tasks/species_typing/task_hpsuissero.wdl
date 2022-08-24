version 1.0

task hpsuissero {
  meta {
    description: "Serotype prediction of Haemophilus parasuis assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/hpsuissero:1.0.1--hdfd78af_0"
    Int? cpu = 4
    String version = "1.0.1"
  }
  command <<<
    # Does not output a version
    echo ~{version} | tee VERSION
    HpsuisSero.sh \
      -i ~{assembly} \
      -o ./ \
      -s ~{samplename} \
      -x fasta \
      -t ~{cpu}
  >>>
  output {
    File hpsuissero_results = "~{samplename}.tsv"
    String hpsuissero_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
