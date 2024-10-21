version 1.0

task pilon {
  input {
    File assembly
    File bam
    File bai
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/pilon:1.24--hdfd78af_0"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # version capture
    pilon --version | cut -d' ' -f3 | tee VERSION

    # run pilon
    if pilon \
    --genome ~{assembly} \
    --frags ~{bam} \
    --output ~{samplename} \
    --outdir pilon \
    --changes --vcf; then
      touch WARNING
    else
      tee "Pilon failed to run for ~{samplename}" > WARNING
      exit 1
    fi

  >>>
  output {
    File? assembly_fasta = "pilon/~{samplename}.fasta"
    File? changes = "pilon/~{samplename}.changes"
    File? vcf = "pilon/~{samplename}.vcf"
    String pilon_version = read_string("VERSION")
    String pilon_docker = "~{docker}"
    String pilon_warning = read_string("WARNING")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}