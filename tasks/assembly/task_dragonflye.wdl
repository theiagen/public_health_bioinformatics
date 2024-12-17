version 1.0

task dragonflye {
  input {
    File read1 # intended for ONT data only
    String samplename
    String? assembler # default is flye
    String? assembler_options # default ''
    String? genome_length # default autodetect
    File? illumina_read1 # illumina read1 fastq to use for polishing
    File? illumina_read2 # illumina read2 fastq to use for polishing 
    Boolean use_pilon_illumina_polisher = false # use polypolish by default
    Int illumina_polishing_rounds = 1 # default 1
    Int polishing_rounds = 1 # default 1
    Boolean use_racon = false # use medaka polishing by default
    String medaka_model = "r941_min_hac_g507"
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/dragonflye:1.0.14--hdfd78af_0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 32
  }
  command <<<
    # get version information
    dragonflye --version | tee VERSION

    # determine polishing parameter
    if [ "~{use_racon}" = true ]; then 
      # --racon indicates number of polishing rounds to conduct with Racon
      POLISHER="--racon ~{polishing_rounds}"
    else
      # --medaka indicates number of polishing rounds to conduct with medaka using medaka_model
      POLISHER="--model ~{medaka_model} --medaka ~{polishing_rounds}"
    fi

    # determine if illumina polishing will be performed based on presence of two files
    if [ -f "~{illumina_read1}" ] && [ -f "~{illumina_read2}" ]; then
      echo "Illumina reads provided; will perform Illumina polishing"
      if [ "~{use_pilon_illumina_polisher}" = true ]; then
        # use --pilon instead of polypolish as default
        ILLUMINA_POLISHER="--R1 ~{illumina_read1} --R2 ~{illumina_read2} --pilon ~{illumina_polishing_rounds}"
      else
        # use --polypolish by default
        ILLUMINA_POLISHER="--R1 ~{illumina_read1} --R2 ~{illumina_read2} --polypolish ~{illumina_polishing_rounds}"
      fi
    else
      echo "Illumina reads were not provided."
      ILLUMINA_POLISHER=""
    fi

    # reduce requested memory to prevent out of memory errors
    mem=$(( ~{memory}-1 ))

    # run dragonflye
    # --reads for input nanopore fastq
    # --depth 0 disables sub-sampling of reads (performed with rasusa already)
    # --outdir indicates output directory
    # --gsize indicates genome size; if this is not set it will autodetect size
    # --cpus number of cpus to use
    # --ram try to keep RAM usage below this number
    # --nofilter disables read length filtering (performed with nanoq already)
    # --assembler has three options: raven, miniasm, flye (default: flye)
    # --opts enables extra assembler options in quotes
    # see above for polisher input explanation
    # see above for illumina polisher input explanation
    dragonflye \
      --reads ~{read1} \
      --depth 0 \
      --outdir dragonflye \
      ~{'--gsize ' + genome_length} \
      --cpus ~{cpu} \
      --ram ${mem} \
      --nofilter \
      ~{'--assembler ' + assembler} \
      ~{'--opts "' + assembler_options + '"'} \
      ${POLISHER} \
      ${ILLUMINA_POLISHER}

    # rename final output file to have .fasta ending instead of .fa
    mv dragonflye/contigs.fa ~{samplename}.fasta
    mv dragonflye/*.gfa ~{samplename}_contigs.gfa
  >>>
  output {
    File assembly_fasta = "~{samplename}.fasta"
    File contigs_gfa = "~{samplename}_contigs.gfa"
    String dragonflye_version = read_string("VERSION")
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