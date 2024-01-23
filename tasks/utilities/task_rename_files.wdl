version 1.0

task rename_PE_files {
  meta {
    description: "Rename a set of Paired-End FASTQ files"
  }
  input {
    File read1
    File read2
    String new_filename
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 2
  }
  command <<<
    # rename forward read
    # check if compressed, if not rename and compress
    if [[ "~{read1}" == *.gz ]]; then
      cp ~{read1} ~{new_filename}_R1.fastq.gz
    else
      cp ~{read1} ~{new_filename}_R1.fastq
      gzip ~{new_filename}_R1.fastq
    fi

    # check if reverse read exists
    if [[ "~{read2}" == *.gz ]]; then
      cp ~{read2} ~{new_filename}_R2.fastq.gz
    else
      cp ~{read2} ~{new_filename}_R2.fastq
      gzip ~{new_filename}_R2.fastq
    fi
  >>>
  output {
    File read1_renamed = "~{new_filename}_R1.fastq.gz"
    File read2_renamed = "~{new_filename}_R2.fastq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task rename_SE_files {
  meta {
    description: "Rename a set of Single-End FASTQ files"
  }
  input {
    File read1
    String new_filename
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 2
  }
  command <<<
    # rename forward read
    # check if compressed, if not rename and compress
    if [[ "~{read1}" == *.gz ]]; then
      cp ~{read1} ~{new_filename}.fastq.gz
    else
      cp ~{read1} ~{new_filename}.fastq
      gzip ~{new_filename}.fastq
    fi
  >>>
  output {
    File read1_renamed = "~{new_filename}.fastq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

