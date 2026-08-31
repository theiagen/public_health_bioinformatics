version 1.0

task allele_calling {
  meta {
    description: "PulseNet 2.0 Allele Calling algorithm"
  }
  input {
    String samplename
    File assembly

    String scheme
    File blast_db
    String loci_path
    Float similarity_threshold
    String qc_genus
    String qc_species

    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/allele_calling:1.0.0-dev"
    Int memory = 32
  }
  command <<<
    # save original directory path for compatibility with BOTH miniwdl and cromwell
    ORIGINAL_DIR=$(pwd)

    # tool fails if not run in /app directory
    cd /app
    tar -xzf ~{blast_db}

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

    # run allele calling algorithm
    ngs-run AlleleCalling \
      --sample-id ~{samplename} \
      --publish-dir . \
      --n-threads ~{cpu} \
      --assembly assembly.fasta.gz \
      --blast-kb.path tests/files/blast_kb \
      --blast-kb.similarity ~{similarity_threshold} \
      --blast-kb.db ~{scheme} \
      --blast-kb.loci ~{loci_path} \
      --qc-kb.path tests/files/qc_kb \
      --organism.genus ~{qc_genus} \
      ~{qc_species}

    # move all results to the original directory for miniwdl/cromwell
    # and then change to that directory for output parsing
    cp *.gz *.json *.bam ${ORIGINAL_DIR}
    cd $ORIGINAL_DIR

    # unzip file for later format conversion
    gunzip calls_core_standard.csv.gz

    # print the outputs json to stdout for debugging purposes
    echo "outputs.json:"
    cat outputs.json

    # parse the results file `outputs.json` for the QC metrics
    # also convert core calls into JSON
    python3 <<CODE
    import json

    # extract the "QC" section from the outputs.json file
    with open("outputs.json") as infile:
      data = json.load(infile)["qc"]

    qc_result = data["result"]
    metrics = data["metrics"]

    # write the qc pass/fail result to a file
    with open("ALLELE_CALLING_RESULT", "w") as outfile:
      outfile.write(str(qc_result))

    # get various qc metrics and write to file
    for field in ("coreCount", "corePercentage", "accessoryCount", "accessoryPercentage", "totalLociCount"):
      with open(field.upper(), "w") as outfile:
        outfile.write(str(metrics[field]))

    ###### convert core allele calls into JSON format ######
    import csv

    with open("calls_core_standard.csv", newline="") as infile:
      reader = csv.DictReader(infile)
      row = next(reader)

    values = {}
    for allele, hash in row.items():
      if allele == "Sample_ID":
        continue
      if hash == "?": # replace ? characters with ""
        values[allele] = ""
      else:
        values[allele] = int(hash)

    with open("calls_core.json", "w") as outfile:
      json.dump({"sample_id": row["Sample_ID"], "values": values}, outfile)

    CODE

    # rezip file
    gzip calls_core.json
  >>>
  output {
    String allele_calling_scheme = scheme
    String allele_calling_result = read_string("ALLELE_CALLING_RESULT")
    File allele_calling_standard_json = "calls_standard.json.gz"
    File allele_calling_core_json = "calls_core.json.gz"
    File allele_calling_bam = "allele_calls.bam"
    # QC metrics
    Int allele_calling_core_count = read_int("CORECOUNT")
    Float allele_calling_core_percentage = read_float("COREPERCENTAGE")
    Int allele_calling_accessory_count = read_int("ACCESSORYCOUNT")
    Float allele_calling_accessory_percentage = read_float("ACCESSORYPERCENTAGE")
    Int allele_calling_total_loci_count = read_int("TOTALLOCICOUNT")
    # Versioning
    String allele_calling_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
