version 1.0

task sra_lite_autodetect {
  meta {
    description: "Detect if SRA-Lite FASTQ file"
  }
  input {
    File read1
    File? read2
    String docker = "us-docker.pkg.dev/general-theiagen/quay/ubuntu:latest"
    Int disk_size = 100
    Int cpu = 1
    Int memory = 2
  }
  command <<<
    # determine if paired-end or not
    if ! [ -z ~{read2} ]; then
      echo "Reads are paired..."
    fi

    # detect if file is compressed and define the appropriate command
    if [[ ~{read1} == *.gz ]]; then
      echo "Reads are compressed..."
      command="zcat"
    else
      echo "Reads are not compressed..."
      command="cat"
    fi

    # check if the first quality control string is set to SRA-Lite
    # SRA-Lite filetype has all the quality enconding set to the '?' character
    # corresponding to a phred-score of 30
    # awk is checking the 4th line of the file and if it starts with and contains only '?' characters
    #   NR==4 on the fourth line,
    #   if $0 ~ /.../ this will evaluate to true only if the entire line matches this regular expression, which is:
    #   ^ at the start of the line
    #   [?] look for this character
    #   +$ and only this character until the end of the line
    $command ~{read1} | head -n 4 | awk 'NR==4 {if ($0 ~ /^[?]+$/) {print "SRA-Lite FASTQ detected"} else {print ""}}' > WARNING
  >>>
  output {
    String warning = read_string("WARNING")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}