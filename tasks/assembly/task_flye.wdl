version 1.0

task flye_consensus {
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

    Int polishing_iterations = 1
    Int minimum_overlap = 1
    
    Float? read_error_rate
    Boolean uneven_coverage_mode = false
    Boolean keep_haplotypes = false
    Boolean no_alt_contigs = false
    Boolean scaffold = false
      
    String? additonal_parameters 

    Int cpu = 4
    Int disk_size = 100
    String docker
    Int memory = 32
  }
  command <<<
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
      --iterations ~{polishing_iterations} \
      --min-overlap ~{minimum_overlap} \
      ~{if defined(asm_coverage) then "--genome-size " + genome_length else ""} \
      ~{"--asm-coverage " + asm_coverage} \
      ~{"--read-error " + read_error_rate} \
      ~{true="--meta" false="" uneven_coverage_mode} \
      ~{true="--keep-haplotypes" false="" keep_haplotypes} \
      ~{true="--no-alt-contigs" false="" no_alt_contigs} \
      ~{true="--scaffold" false="" scaffold} \
      ~{"--extra-params " + additonal_parameters} \
      --threads ~{cpu} \
      --out-dir .

  >>>
  output {
    File assembly = "~{samplename}.assembly.fasta"
    String flye_version = read_string("VERSION")
    String flye_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}