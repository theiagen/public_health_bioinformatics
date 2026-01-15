version 1.0

task stats_n_coverage {
  input {
    File bamfile # aligned reads
    String samplename
    File read1
    File? read2
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
  }
  command <<<
    date | tee DATE
    samtools --version | head -n1 | tee VERSION

    samtools stats ~{bamfile} > ~{samplename}.stats.txt
    samtools coverage ~{bamfile} -m -o ~{samplename}.cov.hist
    samtools coverage ~{bamfile} -o ~{samplename}.cov.txt
    samtools flagstat ~{bamfile} > ~{samplename}.flagstat.txt

     # Extracting coverage, depth, meanbaseq, and meanmapq
    coverage=$(cut -f 6 ~{samplename}.cov.txt | tail -n 1)
    depth=$(cut -f 7 ~{samplename}.cov.txt | tail -n 1)
    meanbaseq=$(cut -f 8 ~{samplename}.cov.txt | tail -n 1)
    meanmapq=$(cut -f 9 ~{samplename}.cov.txt | tail -n 1)

    if [ -z "$coverage" ] ; then coverage="0" ; fi
    if [ -z "$depth" ] ; then depth="0" ; fi
    if [ -z "$meanbaseq" ] ; then meanbaseq="0" ; fi
    if [ -z "$meanmapq" ] ; then meanmapq="0" ; fi

    echo $coverage | tee COVERAGE
    echo $depth | tee DEPTH
    echo $meanbaseq | tee MEANBASEQ
    echo $meanmapq | tee MEANMAPQ

    # parse inputted reads for total read count
    read1_count=$(samtools view -c ~{read1})
    if [ ~{if defined(read2) then "true" else "false"} == "true" ]; then
      read2_count=$(samtools view -c ~{read2})
      total_reads=$(echo $(($read1_count + $read2_count)))
    else
      total_reads=$read1_count
    fi
    # exclude supplementary, unmapped, and secondary alignments from the mapped count
    mapped_reads=$(samtools view -c -F 0x904 ~{bamfile})

    # Check for empty values and set defaults to avoid errors
    if [ -z "$total_reads" ]; then total_reads="1"; fi  # Avoid division by zero
    if [ -z "$mapped_reads" ]; then mapped_reads="0"; fi

    # Calculate the percentage of mapped reads
    percentage_mapped_reads=$(awk "BEGIN {printf \"%.2f\", ($mapped_reads / $total_reads) * 100}")

    # If the percentage calculation fails, default to 0.0
    if [ -z "$percentage_mapped_reads" ]; then percentage_mapped_reads="0.0"; fi

    # Output the result
    echo $percentage_mapped_reads | tee PERCENTAGE_MAPPED_READS

    #output all metrics in one txt file
    # Output header row (for CSV)
    echo -e "Statistic\tValue" > ~{samplename}_metrics.txt

    # Output each statistic as a row
    echo -e "Coverage\t$coverage" >> ~{samplename}_metrics.txt
    echo -e "Depth\t$depth" >> ~{samplename}_metrics.txt
    echo -e "Mean Base Quality\t$meanbaseq" >> ~{samplename}_metrics.txt
    echo -e "Mean Mapping Quality\t$meanmapq" >> ~{samplename}_metrics.txt
    echo -e "Percentage Mapped Reads\t$percentage_mapped_reads" >> ~{samplename}_metrics.txt
  >>>
  output {
    String date = read_string("DATE")
    String samtools_version = read_string("VERSION")
    File stats = "~{samplename}.stats.txt"
    File cov_hist = "~{samplename}.cov.hist"
    File cov_stats = "~{samplename}.cov.txt"
    File flagstat = "~{samplename}.flagstat.txt"
    Float coverage = read_string("COVERAGE")
    Float depth = read_string("DEPTH")
    Float meanbaseq = read_string("MEANBASEQ")
    Float meanmapq = read_string("MEANMAPQ")
    Float percentage_mapped_reads = read_string("PERCENTAGE_MAPPED_READS")
    File metrics_txt = "~{samplename}_metrics.txt"

  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
