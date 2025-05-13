version 1.0

task gamma {
  input {
    File assembly
    File? gamma_db = "gs://theiagen-public-resources-rp/reference_data/databases/gamma/default_ResFinderDB_Combined_05-06-20.fsa"
    Int min_percent_identity = 90
    Int min_length_percent_gammas = 20
    Boolean run_gammas = false
    Boolean output_gff = true
    Boolean output_fasta = true
    Boolean extended_output = false
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/gamma:2.2"
    Int cpu = 2
    Int disk_size = 50
    Int memory = 8
  }
  command <<<
    set -euo pipefail
    echo "~{docker}" | sed 's/.*://' | tee VERSION

    if [ ~{run_gammas} == true ]; then
      GAMMA-S.py ~{assembly} \
        ~{gamma_db} \
        ~{samplename} \
        $(if ~{extended_output}; then echo "-e"; fi) \
        -m ~{min_length_percent_gammas} \
        -i ~{min_percent_identity}
    else
      GAMMA.py ~{assembly} \
        ~{gamma_db} \
        ~{samplename} \
        $(if ~{output_gff}; then echo "--gff"; fi) \
        $(if ~{output_fasta}; then echo "--fasta"; fi) \
        $(if ~{extended_output}; then echo "-e"; fi) \
        -i ~{min_percent_identity} \
        --name 
    fi
  >>>
  output {
    File? gamma_results = "~{samplename}.gamma"
    File? gamma_gff = "~{samplename}.gff"
    File? gamma_fasta = "~{samplename}.fasta"
    String gamma_docker = docker
    String gamma_version = read_string("VERSION")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}