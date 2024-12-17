version 1.0

task flye {
  input {
    File read1
    String samplename

    # data type options; by default, uses --nano-raw
    Boolean ont_corrected = false
    Boolean ont_high_quality = false
    Boolean pacbio_raw = false
    Boolean pacbio_corrected = false
    Boolean pacbio_hifi = false
    
    Int? genome_length # requires `asm_coverage`
    Int? asm_coverage # reduced coverage for initial disjointig assembly

    Int flye_polishing_iterations = 1
    Int? minimum_overlap
    
    Float? read_error_rate
    Boolean uneven_coverage_mode = false
    Boolean keep_haplotypes = false
    Boolean no_alt_contigs = false
    Boolean scaffold = false
      
    String? additional_parameters    # Any extra Flye-specific parameters

    Int cpu = 4
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/flye:2.9.4"
    Int memory = 32
  }
  command <<<
    set -euo pipefail
    flye --version | tee VERSION

    # determine read type
    if ~{ont_corrected}; then
      READ_TYPE="--nano-corr"
    elif ~{ont_high_quality}; then
      READ_TYPE="--nano-hq"
    elif ~{pacbio_raw}; then
      READ_TYPE="--pacbio-raw"
    elif ~{pacbio_corrected}; then
      READ_TYPE="--pacbio-corr"
    elif ~{pacbio_hifi}; then
      READ_TYPE="--pacbio-hifi"
    else
      READ_TYPE="--nano-raw"
    fi

    # genome size parameter requires asm_coverage
    flye \
      ${READ_TYPE} ~{read1} \
      --iterations ~{flye_polishing_iterations} \
      ~{"--min-overlap" + minimum_overlap} \
      ~{if defined(asm_coverage) then "--genome-size " + genome_length else ""} \
      ~{"--asm-coverage " + asm_coverage} \
      ~{"--read-error " + read_error_rate} \
      ~{true="--meta" false="" uneven_coverage_mode} \
      ~{true="--keep-haplotypes" false="" keep_haplotypes} \
      ~{true="--no-alt-contigs" false="" no_alt_contigs} \
      ~{true="--scaffold" false="" scaffold} \
      ~{"--extra-params " + additional_parameters } \
      --threads ~{cpu} \
      --out-dir .

    mv assembly.fasta ~{samplename}.assembly.fasta
    mv assembly_info.txt ~{samplename}.assembly_info.txt
    mv assembly_graph.gfa ~{samplename}.assembly_graph.gfa

  >>>
  output {
    File assembly_fasta = "~{samplename}.assembly.fasta"
    File assembly_graph_gfa = "~{samplename}.assembly_graph.gfa" 
    File assembly_info = "~{samplename}.assembly_info.txt" 
    String flye_version = read_string("VERSION")
    String flye_docker = "~{docker}"
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