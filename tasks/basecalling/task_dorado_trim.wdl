version 1.0

task dorado_trim {
  input {
    Array[File] fastq_files  # FASTQ files from the demultiplexing step
    File custom_primers # Custom primers in FASTA format
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.9.0-cuda12.2.0"
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
  }
  command <<< 
    mkdir -p trimmed_fastqs
    for fq in ~{sep=" " fastq_files}; do
      filename=$(basename "$fq")
      dorado trim "$fq" --primer-sequences "~{custom_primers}" --emit-fastq > trimmed_fastqs/"$filename"
    done
  >>>
  output {
    Array[File] trimmed_fastq_files = glob("trimmed_fastqs/*.fastq.gz")
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
