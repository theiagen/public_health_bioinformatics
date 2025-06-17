version 1.0

task porechop {
  input {
    File read1
    String samplename
    String? trimopts # Optional trimming options
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/porechop:0.2.4"
  }
  command <<< 
    set -euo pipefail

    echo "Porechop version:"
    porechop --version | tee VERSION

    echo "Trimming reads with Porechop..."
    porechop \
      --input ~{read1} \
      --output ~{samplename}.trimmed.fastq \
      --threads ~{cpu} \
      ~{trimopts}

    echo "Compressing trimmed reads..."
    gzip -f ~{samplename}.trimmed.fastq
  >>>
  output {
    File trimmed_reads = "~{samplename}.trimmed.fastq.gz"
    String porechop_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
