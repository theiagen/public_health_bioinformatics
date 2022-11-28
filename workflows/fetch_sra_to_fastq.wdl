version 1.0

workflow fetch_sra_to_fastq {

  input {
    String    SRR
  }

  call prefetch_fastq_dump {
    input:
      sra_id=SRR
  }

  output {
    File    read1   =prefetch_fastq_dump.read1
    File    read2   =prefetch_fastq_dump.read2
  }
}

task prefetch_fastq_dump {

  input {
    String    sra_id
  }

  command {
    prefetch --version | head -1 | tee VERSION
    prefetch ${sra_id}
    fastq-dump --version | head -1 | tee VERSION
    fastq-dump \
    --gzip \
    --split-files \
    ${sra_id}
  }

  output {
    File    read1="${sra_id}_1.fastq.gz"
    File    read2="${sra_id}_2.fastq.gz"
  }

  runtime {
    docker:       "quay.io/staphb/sratoolkit:2.9.2"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  1
  }
}
