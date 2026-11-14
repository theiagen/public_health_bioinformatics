version 1.0

task nanoplot {
  input {
    File read1 # intended for ONT data only
    String samplename
    Int max_length = 100000
    Int? est_genome_length

    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/nanoplot:1.48.0-seqkit-dev"
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    # get version
    NanoPlot --version | tee "VERSION"
    seqkit version | tee SEQKIT_VERSION

    # run nanoplot
    # --prefix for output file tag
    # --threads for number of threads allowed
    # --N50 to display N50 mark in read length histogram
    # --loglength to show logarithmic scaling of lengths
    # --tsv_stats to output the stats file in TSV format
    # --maxlength to hide reads longer than this
    # --no_static to skip static PNG plots, which require a runtime Chrome download and are not task outputs
    NanoPlot \
      --fastq ~{read1} \
      --prefix "~{samplename}_" \
      --threads 4 \
      --N50 \
      --loglength \
      --tsv_stats \
      --maxlength ~{max_length} \
      --no_static

    # only populate read statistics, GC content, and estimated coverage if NanoPlot produced a non-empty stats file
    if [[ -s ~{samplename}_NanoStats.txt ]]; then
      # grep read statistics from tsv stats file
      grep "number_of_reads" ~{samplename}_NanoStats.txt | cut -f 2 | tee NUMBER_OF_READS
      NUM_BASES=$(grep "number_of_bases" ~{samplename}_NanoStats.txt | cut -f 2 | tee NUMBER_OF_BASES)
      grep "median_read_length" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEDIAN_READ_LENGTH
      grep "mean_read_length" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEAN_READ_LENGTH
      grep "read_length_stdev" ~{samplename}_NanoStats.txt | cut -f 2 | tee READ_LENGTH_STDEV
      grep "n50" ~{samplename}_NanoStats.txt | cut -f 2 | tee N50
      grep "mean_qual" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEAN_QUAL
      grep "median_qual" ~{samplename}_NanoStats.txt | cut -f 2 | tee MEDIAN_QUAL

      # capture GC content
      echo "DEBUG: running seqkit stats on $(basename ~{read1})"
      seqkit stats -a -T -j ~{cpu} ~{read1} | tee seqkit_stats.tsv
      # column 18 is GC(%) in seqkit v2.14.0 `stats -a -T` output; revisit if the docker's seqkit version changes
      awk -F'\t' 'NR==2 {print $18}' seqkit_stats.tsv > GC_PERCENT
      echo "DEBUG: GC content in $(basename ~{read1}): $(cat GC_PERCENT)"
      # add GC content to the NanoStats file alongside the other summary statistics
      awk -F'\t' -v OFS='\t' -v gc="$(cat GC_PERCENT)" '{print} $1 == "median_qual" {print "gc_percent", gc}' \
        ~{samplename}_NanoStats.txt > tmp.txt
      mv tmp.txt ~{samplename}_NanoStats.txt

      if [ ~{if defined(est_genome_length)  then "true" else "false" } == true ]; then
        # estimate coverage
        # using math: C = N / G where N is number of bases, and G is estimated genome size
        python3 -c "print(round(${NUM_BASES} / ~{est_genome_length}, 2))" | tee EST_COVERAGE
      fi
    fi

  >>>
  output {
    File? nanoplot_html = "~{samplename}_NanoPlot-report.html"
    File? nanoplot_tsv = "~{samplename}_NanoStats.txt"
    Int? num_reads = read_int("NUMBER_OF_READS")
    Float? median_readlength = read_float("MEDIAN_READ_LENGTH")
    Float? mean_readlength = read_float("MEAN_READ_LENGTH")
    Float? stdev_readlength = read_float("READ_LENGTH_STDEV")
    Float? n50 = read_float("N50")
    Float? mean_q = read_float("MEAN_QUAL")
    Float? median_q = read_float("MEDIAN_QUAL")
    Float? gc_percent = read_float("GC_PERCENT")
    Float? est_coverage = read_float("EST_COVERAGE")
    String nanoplot_version = read_string("VERSION")
    String seqkit_version = read_string("SEQKIT_VERSION")
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
