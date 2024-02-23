version 1.0

import "../../../tasks/utilities/task_sra_fetch.wdl" as sra_fetch

workflow fetch_sra_to_fastq {
  input {
    String sra_accession
    String? docker
    Int? disk_size
    Int? memory
    Int? cpu
    String? fastq_dl_opts
  }
  call sra_fetch.fastq_dl_sra {
    input:
      sra_accession = sra_accession,
      docker = docker,
      disk_size = disk_size,
      cpu = cpu,
      memory = memory,
      fastq_dl_opts = fastq_dl_opts
  }
  output {
    File read1 = fastq_dl_sra.read1
    File? read2 = fastq_dl_sra.read2
    File fastq_dl_fastq_metadata = fastq_dl_sra.fastq_metadata
    String fastq_dl_version = fastq_dl_sra.fastq_dl_version
    String fastq_dl_docker = fastq_dl_sra.fastq_dl_docker
    String fastq_dl_date = fastq_dl_sra.fastq_dl_date
  }
}