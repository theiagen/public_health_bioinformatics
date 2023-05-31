version 1.0

task samtools_mpileup {
  meta {
    description: "Produce a pileup file with samtools mpileup"
  }
  input {
    File bam
    File bai
    File reference
    String samplename
    String docker = "staphb/samtools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    samtools mpileup -f ~{reference} ~{bam} > "~{samplename}".pileup
  >>>
  output {
    File samtools_pileup = "~{samplename}.pileup"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}