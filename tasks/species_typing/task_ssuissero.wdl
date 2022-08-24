version 1.0

task ssuissero {
  meta {
    description: "Serotype prediction of Streptococcus suis assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/ssuissero:1.0.1--hdfd78af_0"
    Int? cpu = 4
    String version = "1.0.1"
  }
  command <<<
    # Does not output a version
    echo ~{version} | tee VERSION
    SsuisSero.sh \
      -i ~{assembly} \
      -o ./ \
      -s ~{samplename} \
      -x fasta \\
      -t ~{cpu}
  >>>
  output {
    File ssuissero_results = "~{samplename}.tsv"
    String ssuissero_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
