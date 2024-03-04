version 1.0

task hpsuissero {
  meta {
    description: "Serotype prediction of Haemophilus parasuis assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/hpsuissero:1.0.1--hdfd78af_0"
    Int cpu = 4
    String version = "1.0.1"
    Int disk_size = 50
    Int memory = 8
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
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
