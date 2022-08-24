version 1.0

task bbduk_pe {
  input {
    File read1_trimmed
    File read2_trimmed
    String samplename
    Int mem_size_gb=8
    String docker = "quay.io/staphb/bbtools:38.76"
  }
  command <<<
    # date and version control
    date | tee DATE

    repair.sh in1=~{read1_trimmed} in2=~{read2_trimmed} out1=~{samplename}.paired_1.fastq.gz out2=~{samplename}.paired_2.fastq.gz

    bbduk.sh in1=~{samplename}.paired_1.fastq.gz in2=~{samplename}.paired_2.fastq.gz out1=~{samplename}.rmadpt_1.fastq.gz out2=~{samplename}.rmadpt_2.fastq.gz ref=/bbmap/resources/adapters.fa stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

    bbduk.sh in1=~{samplename}.rmadpt_1.fastq.gz in2=~{samplename}.rmadpt_2.fastq.gz out1=~{samplename}_1.clean.fastq.gz out2=~{samplename}_2.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=~{samplename}.phix.stats.txt

  >>>
  output {
    File read1_clean = "~{samplename}_1.clean.fastq.gz"
    File read2_clean = "~{samplename}_2.clean.fastq.gz"
    File adapter_stats = "~{samplename}.adapters.stats.txt"
    File phiX_stats = "~{samplename}.phix.stats.txt"
    String bbduk_docker = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{mem_size_gb} GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}

task bbduk_se {
  input {
    File read1_trimmed
    String samplename
    Int mem_size_gb=8
    String docker="quay.io/staphb/bbtools:38.76"
  }
  command <<<
    # date and version control
    date | tee DATE

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}.rmadpt_1.fastq.gz ref=/bbmap/resources/adapters.fa stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}_1.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=~{samplename}.phix.stats.txt
  >>>
  output {
    File read1_clean = "~{samplename}_1.clean.fastq.gz"
    File adapter_stats = "~{samplename}.adapters.stats.txt"
    File phiX_stats = "~{samplename}.phix.stats.txt"
    String bbduk_docker = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{mem_size_gb} GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}
