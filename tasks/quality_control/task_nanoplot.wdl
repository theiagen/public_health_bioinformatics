version 1.0

task nanoplot {
  input {
    File read1 # intended for ONT data only
    String samplename
    Int max_length = 100000
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0"
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    # get version
    NanoPlot --version | tee "VERSION"

    # run nanoplot
    # --prefix for output file tag
    # --threads for number of threads allowed
    # --N50 to display N50 mark in read length histogram
    # --loglength to show logarithmic scaling of lengths
    # --tsv_stats to output the stats file in TSV format
    # --maxlength to hide reads longer than this
    NanoPlot \
      --fastq ~{read1} \
      --prefix "~{samplename}_" \
      --threads 4 \
      --N50 \
      --loglength \
      --tsv_stats \
      --maxlength ~{max_length}
   
  >>>
  output {
    File nanoplot_html = "~{samplename}_NanoPlot-report.html"
    File nanoplot_tsv = "~{samplename}_NanoStats.txt"
    String nanoplot_version = read_string("VERSION")
    String nanoplot_docker = docker
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}