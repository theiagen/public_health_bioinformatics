version 1.0

task nanoplot {
  input {
    File read1 # intended for ONT data only
    String samplename
    Int max_length = 100000
    Int disk_size = 100
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
    String nanoplot_version = read_string("VERSION")
  }
  runtime {
    docker: "quay.io/staphb/nanoplot:1.40.0"
    memory: "16 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}