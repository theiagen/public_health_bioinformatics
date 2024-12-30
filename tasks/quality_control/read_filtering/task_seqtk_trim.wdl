version 1.0

task seqtk_trim {
  input {
    File fastq_file
    String sample_name
    Int left_trim = 25    # Default: trim 25 bp from the left
    Int right_trim = 25   # Default: trim 25 bp from the right
    Int min_length = 50   # Default: minimum read length
    Int memory = 4
    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqtk:1.4"
  }
  command <<< 
    # Trim and filter reads using seqtk
    seqtk trimfq -b ~{left_trim} -e ~{right_trim} ~{fastq_file} > trimmed.fq
    seqtk seq -L ~{min_length} trimmed.fq > ~{sample_name}_filtered.fq
  >>>
  output {
    File trimmed_fastq = "~{sample_name}_filtered.fq"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
