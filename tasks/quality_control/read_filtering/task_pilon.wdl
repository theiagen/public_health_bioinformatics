version 1.0

task pilon {
  input {
    File assembly
    File bam
    File bai
    String samplename
    Int? min_mapping_quality
    Int? min_base_quality
    Float? min_depth
    String? fix
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/pilon:1.24--hdfd78af_0"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # version capture
    pilon --version | cut -d' ' -f3 | tee VERSION

    # run pilon
    pilon \
    --genome ~{assembly} \
    --frags ~{bam} \
    --output ~{samplename} \
    --outdir pilon \
    --changes --vcf \
    ~{'--minmq ' + min_mapping_quality} \
    ~{'--minqual ' + min_base_quality} \
    ~{'--mindepth ' + min_depth} \
    ~{'--fix ' + fix}

  >>>
  output {
    File assembly_fasta = "pilon/~{samplename}.fasta"
    File changes = "pilon/~{samplename}.changes"
    File vcf = "pilon/~{samplename}.vcf"
    String pilon_version = read_string("VERSION")
    String pilon_docker = "~{docker}"
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