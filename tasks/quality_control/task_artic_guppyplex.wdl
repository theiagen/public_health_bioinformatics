version 1.0

task read_filtering {
  input {
    File demultiplexed_reads
    String samplename
    String run_prefix = "artic_ncov2019"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/artic-ncov2019:1.3.0-medaka-1.4.3"
    Int min_length = 400
    Int max_length = 700
    Int cpu = 8
    Int disk_size = 100
    Int memory = 16
  }
  command <<<
    # date and version control
    mkdir ~{samplename}
    cp ~{demultiplexed_reads} ~{samplename}/
    echo "DIRNAME: $(dirname)"

    # run artic guppyplex 
    artic guppyplex --min-length ~{min_length} --max-length ~{max_length} --directory ~{samplename} --prefix ~{run_prefix}
  >>>
  output {
    File filtered_reads = "~{run_prefix}_~{samplename}.fastq"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}