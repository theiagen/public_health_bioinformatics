version 1.0

task ncbi_scrub_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int memory = 8 
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE

    # detect compression and define cat command
    if [[ "~{read1}" == *.gz ]]
    then
      cat_command="zcat"
    else
      cat_command="cat"
    fi

    # if compressed, unzip read files as scrub tool does not take in .gz fastq files, and interleave them
    # paste command takes 4 lines at a time and merges them into a single line with tabs
    # tr substitutes the tab separators from paste into new lines, effectively interleaving the reads and keeping the FASTQ format
    paste <($cat_command ~{read1} | paste - - - -) <($cat_command ~{read2} | paste - - - -) | tr '\t' '\n' > interleaved.fastq

    # dehost reads
    # -x Remove spots instead of default 'N' replacement.
    # -s ; Input is (collated) interleaved paired-end(read) file AND you wish both reads masked or removed.
    # |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED ; capture the number of spots removed by fetching the first item on the last line of stderr
    /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -s -i interleaved.fastq |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED

    # split interleaved reads and compress files
    # paste takes 8 at a time and merges them into a single line with tabs
    # grouping each pair of reads into one line
    paste - - - - - - - - < interleaved.fastq.clean \
      | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > ~{samplename}_R1_dehosted.fastq.gz) \
      | cut -f 5-8 | tr '\t' '\n' | gzip > ~{samplename}_R2_dehosted.fastq.gz
      
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    File read2_dehosted = "~{samplename}_R2_dehosted.fastq.gz"
    Int human_spots_removed = read_int("SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
  }
  runtime {
      docker: "~{docker}"
      memory: memory + " GB"
      cpu: cpu
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}

task ncbi_scrub_se {
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1"
    Int disk_size = 100
    Int memory = 8 
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE

    # unzip fwd file as scrub tool does not take in .gz fastq files
    if [[ "~{read1}" == *.gz ]]
    then
      gunzip -c ~{read1} > r1.fastq
      read1_unzip=r1.fastq
    else
      read1_unzip=~{read1}
    fi

    # dehost reads
    # -x Remove spots instead of default 'N' replacement.
    # |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED ; capture the number of spots removed by fetching the first item on the last line of stderr
    /opt/scrubber/scripts/scrub.sh -p ~{cpu} -x -i ${read1_unzip} |& tail -n1 | awk -F" " '{print $1}' > SPOTS_REMOVED

    # gzip dehosted reads
    gzip ${read1_unzip}.clean -c > ~{samplename}_R1_dehosted.fastq.gz
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    Int human_spots_removed = read_int("SPOTS_REMOVED")
    String ncbi_scrub_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}