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
    # tool fails if not run in /app directory
    cd /app

    tar -xzf ~{blast_db}

    # parse genus and species inputs
    GENUS=$(echo ~{organism} | cut -f1 -d' ')
    SPECIES=$(echo ~{organism} | cut -f2 -d' ')
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

    # check if fasta is gzipped or not; then copy to standard name
    if gzip -t ~{assembly}; then
      # unzip assembly to enable stripping headers of extra content
      gunzip ~{assembly} -c > assembly.fasta
    else
      cp ~{assembly} assembly.fasta
    fi

    # strip headers of extra content and rezip
    sed -i '/^>/s/[[:space:]].*$//' assembly.fasta
    gzip assembly.fasta

    echo $SPECIES $GENUS

    # run allele calling algorithm
    if ngs-run AlleleCalling \
        --sample-id ~{samplename} \
        --publish-dir . \
        --n-threads ~{cpu} \
        --assembly assembly.fasta.gz \
        --blast-kb.path tests/files/blast_kb \
        --blast-kb.similarity ~{similarity_threshold} \
        --blast-kb.db ${SCHEME} \
        --blast-kb.loci ${LOCI_LIST} \
        --qc-kb.path tests/files/qc_kb \
        --organism.genus ${SCHEME}; then
      echo "PASS" > ALLELE_CALLER_RESULT

      cp * /mnt/miniwdl_task_container/work/
      cd /mnt/miniwdl_task_container/work/
    else
      echo "FAIL" > ALLELE_CALLER_RESULT
    fi


    cat outputs.json

    # parse outputs.json for the QC metrics
    python3 <<CODE
    import json

    with open("outputs.json") as infile:
      metrics = json.load(infile)["qc"]["metrics"]

    for field in ("coreCount", "corePercentage", "accessoryCount", "accessoryPercentage", "totalLociCount"):
      with open(field.upper(), "w") as outfile:
        outfile.write(str(metrics[field]))

    CODE

  >>>
  output {
    String allele_caller_scheme = read_string("ALLELE_CALLER_SCHEME")
    String allele_caller_result = read_string("ALLELE_CALLER_RESULT")
    File allele_caller_wgmlst_json = "calls_standard.json.gz"
    File allele_caller_cgmlst_json = "calls_core_standard.csv.gz"
    File allele_caller_detailed_json = "allele_calls.json.gz"
    # QC metrics
    Int allele_caller_core_count = read_int("CORECOUNT")
    Float allele_caller_core_percentage = read_float("COREPERCENTAGE")
    Int allele_caller_accessory_count = read_int("ACCESSORYCOUNT")
    Float allele_caller_accessory_percentage = read_float("ACCESSORYPERCENTAGE")
    Int allele_caller_total_loci_count = read_int("TOTALLOCICOUNT")
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
