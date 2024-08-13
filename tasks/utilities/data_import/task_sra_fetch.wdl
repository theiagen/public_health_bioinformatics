version 1.0


task fastq_dl_sra {
  input {
    String sra_accession
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
    # default set to force the use of SRA instead of ENA due to SRA Lite FASTQ file format issues
    String fastq_dl_opts = "--provider sra"
  }
  meta {
    # so that call caching is always turned off
    volatile: true
  }
  command <<<
    # capture version
    fastq-dl --version | tee VERSION

    # capture date in UTC timezone
    date -u | tee DATE

    # download fastq files
    fastq-dl \
      --verbose \
      -a ~{sra_accession} \
      --cpus ~{cpu} \
      --prefix ~{sra_accession} \
      ~{fastq_dl_opts}

    # tag single-end reads with _1
    if [ -f "~{sra_accession}.fastq.gz" ] && [ ! -f "~{sra_accession}_1.fastq.gz" ]; then
      mv "~{sra_accession}.fastq.gz" "~{sra_accession}_1.fastq.gz"
    fi
    
    # check if the first quality control string is set to SRA-Lite
    # SRA-Lite filetype has all the quality enconding set to the '?' or '$' character
    # corresponding to a phred-score of 30 or 3
    # awk is checking the 4th line of the file and if it starts with and contains only '?' or '$' characters
    #   NR==4 on the fourth line,
    #   if $0 ~ /.../ this will evaluate to true only if the entire line matches this regular expression, which is:
    #   ^ at the start of the line
    #   [?] look for this character
    #   +$ and only this character until the end of the line
    zcat "~{sra_accession}_1.fastq.gz" | head -n 4 | awk 'NR==4 {if ($0 ~ /^[?$]+$/) {print "Potential SRA-Lite FASTQ detected"} else {print ""}}' > WARNING
    
  >>>
  output {
    File read1 = "~{sra_accession}_1.fastq.gz"
    File? read2 = "~{sra_accession}_2.fastq.gz"
    File fastq_metadata = "~{sra_accession}-run-info.tsv"
    String fastq_dl_version = read_string("VERSION")
    String fastq_dl_docker = docker
    String fastq_dl_date = read_string("DATE")
    String warning = read_string("WARNING")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible:  1
    maxRetries: 3
  }
}
