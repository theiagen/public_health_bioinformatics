version 1.0

task filter_contigs {
  input {
    String samplename
    File assembly_fasta
    Int min_length = 1000
    Float min_coverage = 2.0
    Boolean skip_length_filter = false
    Boolean skip_coverage_filter = false
    Boolean skip_homopolymer_filter = false
    Int disk_size = 50
    Int memory = 8
    Int cpu = 1
    # Found in theiagen-docker-builds/assembly-filter
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/shovilter:0.1"
  }
  command <<< 
    set -euo pipefail
    
    echo "Filtering contigs from ~{assembly_fasta}" >&2

    python /scripts/assembly-shovilter.py \
      -i ~{assembly_fasta} \
      -o ~{samplename}_filtered_contigs.fasta \
      -m ~{samplename}_filtering_metrics.txt \
      --minlen ~{min_length} \
      --mincov ~{min_coverage} \
      ~{true="--skip-length-filter" false="" skip_length_filter} \
      ~{true="--skip-coverage-filter" false="" skip_coverage_filter} \
      ~{true="--skip-homopolymer-filter" false="" skip_homopolymer_filter}

  >>>
  output {
    File filtered_fasta = "~{samplename}_filtered_contigs.fasta"
    File assembly_filtering_metrics = "~{samplename}_filtering_metrics.txt"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 3
    preemptible: 0
  }
}
