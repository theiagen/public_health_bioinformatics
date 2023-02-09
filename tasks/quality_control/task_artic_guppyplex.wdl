version 1.0

task read_filtering { # in use
  input {
    File demultiplexed_reads
    String samplename
    String run_prefix = "artic_ncov2019"
    Int min_length = 400
    Int max_length = 700
    Int cpu = 8
    Int disk_size = 100
  }
  command <<<
    # date and version control
    mkdir ~{samplename}
    cp ~{demultiplexed_reads} ~{samplename}/
    echo "DIRNAME: $(dirname)"
    artic guppyplex --min-length ~{min_length} --max-length ~{max_length} --directory ~{samplename} --prefix ~{run_prefix}

  >>>
  output {
    File filtered_reads = "~{run_prefix}_~{samplename}.fastq"
  }
  runtime {
    docker: "quay.io/staphb/artic-ncov2019:1.3.0-medaka-1.4.3"
    memory: "16 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}