version 1.0

task dragonflye {
  input {
    File reads
    String samplename
    String genome_size
    String docker = "quay.io/biocontainers/dragonflye:1.0.14--hdfd78af_0"
    Int disk_size = 100
    Int cpu = 4
  }
  command <<<

    dragonflye --version | cut -d ' ' -f2 | tee VERSION

    dragonflye \
      --reads "~{reads}" \
      --outdir dragonflye \
      --depth 0 \
      --gsize "~{genome_size}" \
      --prefix "~{samplename}" \
      --cpus "~{cpu}" \
      --ram 8 \
      --nofilter

  >>>
  output {
    File assembly_fasta = "dragonflye/~{samplename}.fa"
    File assembly_gfa = "dragonflye/~{samplename}.gfa"
    String dragonflye_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}