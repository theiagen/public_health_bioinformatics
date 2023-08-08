version 1.0

task hostile_pe {
  input {
    File read1
    File read2
    String aligner = "bowtie2"
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/hostile:0.1.0--pyhdfd78af_0"
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    hostile --version | tee Version

    # run hostile
    hostile clean  \
      --fastq1 ~{read1} \
      --fastq2 ~{read2} \
      --aligner ~{aligner}

    # rename output reads - to do: fastq.gz or fq.gz termination
    filename_without_extension_1=$(basename "~{read1}" .fastq.gz)
    filename_without_extension_2=$(basename "~{read2}" .fastq.gz)

    mv "${filename_without_extension_1}.clean_1.fastq.gz" "~{samplename}_R1_dehosted.fastq.gz"
    mv "${filename_without_extension_2}.clean_2.fastq.gz" "~{samplename}_R2_dehosted.fastq.gz"
    
  >>>
  output {
    File read1_dehosted = "~{samplename}_R1_dehosted.fastq.gz"
    File read2_dehosted = "~{samplename}_R2_dehosted.fastq.gz"
    String hostile_docker = docker

  }
  runtime {
      docker: "~{docker}"
      memory: "8 GB"
      cpu: 4
      disks:  "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB" # TES
      preemptible: 0
      maxRetries: 3
  }
}
