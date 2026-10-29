version 1.0

task nanoplot {
  input {
    File read1 # intended for ONT data only
    String samplename
    Int max_length = 100000
    Int? est_genome_length

    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/nanoplot:1.40.0"
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    # get version
    NanoPlot --version | tee "VERSION"

    if [[ $(zcat -f ~{read1} | head -n 4 | wc -l) -lt 4 ]]; then
      echo "FAIL: no reads in input" | tee STATUS
      for metric in NUMBER_OF_READS NUMBER_OF_BASES MEDIAN_READ_LENGTH MEAN_READ_LENGTH READ_LENGTH_STDEV N50 MEAN_QUAL MEDIAN_QUAL EST_COVERAGE; do
        echo 0 > "$metric"
      done
      echo -e "Metrics\tdataset\nnumber_of_reads\t0\nnumber_of_bases\t0" > ~{samplename}_NanoStats.txt
      echo "<html><body>No reads in input; NanoPlot was not run.</body></html>" > ~{samplename}_NanoPlot-report.html
      exit 0
    fi
    echo "PASS" > STATUS

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

    # grep read statistics from tsv stats file
    grep "number_of_reads" ~{samplename}_NanoStats.txt | cut -f 2 | tee NUMBER_OF_READS
    NUM_BASES=$(grep "number_of_bases" ~{samplename}_NanoStats.txt | cut -f 2 | tee NUMBER_OF_BASES)
    grep "median_read_length" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEDIAN_READ_LENGTH
    grep "mean_read_length" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEAN_READ_LENGTH
    grep "read_length_stdev" ~{samplename}_NanoStats.txt | cut -f 2 | tee READ_LENGTH_STDEV
    grep "n50" ~{samplename}_NanoStats.txt | cut -f 2 | tee N50
    grep "mean_qual" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEAN_QUAL
    grep "median_qual" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEDIAN_QUAL

    if [ ~{if defined(est_genome_length)  then "true" else "false" } == true ]; then
      # estimate coverage
      # using math: C = N / G where N is number of bases, and G is estimated genome size
      python3 -c "print(round(${NUM_BASES} / ~{est_genome_length}, 2))" | tee EST_COVERAGE
    fi

  >>>
  output {
    File nanoplot_html = "~{samplename}_NanoPlot-report.html"
    File nanoplot_tsv = "~{samplename}_NanoStats.txt"
    Int num_reads = read_int("NUMBER_OF_READS")
    Float median_readlength = read_float("MEDIAN_READ_LENGTH")
    Float mean_readlength = read_float("MEAN_READ_LENGTH")
    Float stdev_readlength = read_float("READ_LENGTH_STDEV")
    Float n50 = read_float("N50")
    Float mean_q = read_float("MEAN_QUAL")
    Float median_q = read_float("MEDIAN_QUAL")
    Float? est_coverage = read_float("EST_COVERAGE")
    String nanoplot_status = read_string("STATUS")
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
