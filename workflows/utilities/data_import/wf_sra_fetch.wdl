version 1.0

workflow fetch_sra_to_fastq {
  input {
    String sra_accession
    String? docker
    Int? disk_size
    Int? memory
    Int? cpus
  }
  call fastq_dl_sra {
    input:
      sra_accession = sra_accession,
      docker = docker,
      disk_size = disk_size,
      cpus = cpus,
      memory = memory
  }
  output {
    File read1 = fastq_dl_sra.read1
    File? read2 = fastq_dl_sra.read2
  }
}

task fastq_dl_sra {
  input {
    String sra_accession
    String docker = "quay.io/biocontainers/fastq-dl:2.0.1--pyhdfd78af_0"
    Int disk_size = 100
    Int cpus = 2
    Int memory = 8
  }
  command <<<
    fastq-dl --version | tee VERSION
    fastq-dl -a ~{sra_accession}

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
    docker: docker
    memory: memory + " GB"
    cpu: cpus
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible:  1
  }
}
