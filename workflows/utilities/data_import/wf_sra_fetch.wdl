version 1.0

import "../../../tasks/utilities/task_sra_fetch.wdl" as sra_fetch

workflow fetch_sra_to_fastq {
  input {
    String sra_accession
    String? docker
    Int? disk_size
    Int? memory
    Int? cpus
    String? fastq_dl_opts
  }
  call sra_fetch.fastq_dl_sra {
    input:
      sra_accession = sra_accession,
      docker = docker,
      disk_size = disk_size,
      cpus = cpus,
      memory = memory,
      fastq_dl_opts = fastq_dl_opts
  }
  output {
    File read1 = fastq_dl_sra.read1
    File? read2 = fastq_dl_sra.read2
  }
}