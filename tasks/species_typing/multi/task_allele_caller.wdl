version 1.0

task allele_caller {
  meta {
    description: "PulseNet 2.0 Allele Calling algorithm"
  }
  input {
    String samplename
    File assembly

    String organism
    String? subspecies

    File blast_db
    Float similarity_threshold

    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/pulsenet2-0/allelecaller:1.0.0"
    Int memory = 32
  }
  command <<<
    SCHEME=$(basename ~{blast_db} .tar.gz | tee ALLELE_CALLER_SCHEME)
    cd /app

    tar -xzf ~{blast_db}

    GENUS=$(echo ${organism} | cut -f1 -d' ')
    SPECIES=$(echo ${organism} | cut -f2 -d' ')
    LOCI_LIST="${SCHEME}/loci.tsv"

    if [[ $GENUS == $SPECIES ]]; then
      # organism was only genus-level resolution; set SPECIES to blank
      SPECIES=""
    elif [[ ${GENUS} == "Vibrio" ]]; then
      # species-specific loci list for vibrio
      # loci list is in a separate location for vibrio and has a different naming format
      # genus is capitalized with `^^`, which is bash 4.0+ parameter expansion
      if [[ $SPECIES == "parahaemolyticus" || $SPECIES == "cholerae" ]]; then
        LOCI_LIST="${SCHEME}/loci/${GENUS^^}_${SPECIES}_loci.tsv"
      else
        LOCI_LIST="${SCHEME}/loci/${GENUS^^}_loci.tsv"
      fi
    fi

    # check if fasta is gzipped or not
    if gzip -t ~{assembly}; then
      cp ~{assembly} assembly.fasta.gz
    else
      cp ~{assembly} assembly.fasta
      gzip assembly.fasta
    fi

    ngs-run AlleleCalling \
      --sample-id ~{samplename} \
      --publish-dir . \
      --n-threads ~{cpu} \
      --assembly assembly.fasta.gz \
      --blast-kb.path tests/files/blast_kb \
      --blast-kb.similarity ~{similarity_threshold} \
      --blast-kb.db ${SCHEME} \
      --blast-kb.loci ${LOCI_LIST} \
      --qc-kb.path tests/files/qc_kb \
      --organism.genus ${GENUS} \
      --organism.species ${SPECIES}

ngs-run AlleleCalling \
  --sample-id campy \
  --publish-dir . \
  --assembly campylobacter_jejuni.fasta.gz \
  --blast-kb.similarity 70 \
  --blast-kb.path tests/files/blast_kb \
  --blast-kb.db CAMPY \
  --blast-kb.loci CAMPY/loci.tsv \
  --qc-kb.path tests/files/qc_kb/ \
  --organism.genus Campylobacter \
  --organism.species jejuni

  echo "test" > ALLELE_CALLER_RESULT


  >>>
  output {
    String allele_caller_scheme = read_string("ALLELE_CALLER_SCHEME")
    String allele_caller_result = read_string("ALLELE_CALLER_RESULT")
    File allele_caller_wgmlst_json = "calls_standard.json.gz"
    File allele_caller_cgmlst_json = "calls_core_standard.json.gz"
    File allele_caller_detailed_json = "allele_calls.json.gz"
    # QC metrics
    # Int allele_caller_core_count = read_int("CORE_COUNT")
    # Float allele_caller_core_percentage = read_float("CORE_PERCENTAGE")
    # Int allele_caller_accessory_count = read_int("ACCESSORY_COUNT")
    # Float allele_caller_accessory_percentage = read_float("ACCESSORY_PERCENTAGE")
    # Int allele_caller_total_loci_count = read_int("TOTAL_LOCI_COUNT")
    # Versioning
    # String allele_caller_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
 #   maxRetries: 3
    preemptible: 1
  }
}
