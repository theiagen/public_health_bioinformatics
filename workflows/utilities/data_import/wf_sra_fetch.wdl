version 1.0

workflow fetch_sra_to_fastq {
  input {
    String sra_accession
  }
  call fastq_dl_sra {
    input:
      sra_accession=sra_accession
  }
  output {
    File read1 = fastq_dl_sra.read1
    File? read2 = fastq_dl_sra.read2
  }
}

task fastq_dl_sra {
  input {
    String sra_accession
  }
  command <<<
    fastq-dl --version | tee VERSION
    fastq-dl ~{sra_accession} SRA

    # tag single-end reads with _1
    if [ -f "~{sra_accession}.fastq.gz" ] && [ ! -f "~{sra_accession}_1.fastq.gz" ]; then
      mv "~{sra_accession}.fastq.gz" "~{sra_accession}_1.fastq.gz"
    fi
  >>>
  output {
    File read1 = "~{sra_accession}_1.fastq.gz"
    File? read2 = "~{sra_accession}_2.fastq.gz"
  }
  runtime {
    docker: "quay.io/biocontainers/fastq-dl:1.1.0--hdfd78af_0"
    memory:"8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible:  1
  }
}
