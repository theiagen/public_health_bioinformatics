version 1.0

task dorado_trim {
  input {
    Array[File] fastq_files  # FASTQ files from the demultiplexing step
    File custom_primers # Custom primers in FASTA format
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.3"
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
  }
  command <<< 
    dorado trim ~{sep=" " fastq_files} --primer-sequences ~{custom_primers} --emit-fastq
  >>>
  output {
    Array[File] trimmed_fastq_files = glob("*.fastq.gz")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
