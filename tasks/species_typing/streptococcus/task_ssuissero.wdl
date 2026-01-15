version 1.0

task ssuissero {
  meta {
    description: "Serotype prediction of Streptococcus suis assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/ssuissero:1.0.1--hdfd78af_0"
    Int cpu = 4
    String version = "1.0.1"
    Int disk_size = 100
    Int memory = 8
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
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
