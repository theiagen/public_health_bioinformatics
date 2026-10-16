version 1.0

task seqsender_build_and_prep {
  input {
    String terra_project_name
    String terra_workspace_name
    String terra_table_name
    File seqsender_config # NCBI password already encrypted with `seqsender.py load_credentials`
    Array[String] sample_ids

    String bioproject
    Array[File]? assembly_fastas # required for GENBANK
    Array[File]? read1_fastqs # required for SRA
    Array[File]? read2_fastqs # optional for SRA (single-end)
    String? sra_transfer_gcp_bucket # NCBI-readable bucket the passing samples' reads are copied to

    String organism
    Boolean submit_to_biosample
    Boolean submit_to_sra
    Boolean submit_to_genbank

    Boolean prep_only

    File? input_table # for command line testing only; used instead of downloading the Terra table

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/seqsender:1.5.2-dev"
    Int cpu = 2
    Int memory = 8
    Int disk_size = 25
  }
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail

    micromamba run -n seqsender \
    seqsender.py version | grep -o "[0-9][0-9.]*[0-9]" | head -n 1 | tee VERSION || true

    # bash arrays; build_and_prep.py checks that the file arrays line up with the sample IDs
    SAMPLE_IDS=(~{sep=' ' sample_ids})
    ASSEMBLY_FASTAS=(~{sep=' ' assembly_fastas})
    READ1_FASTQS=(~{sep=' ' read1_fastqs})
    READ2_FASTQS=(~{sep=' ' read2_fastqs})

    # download the Terra {terra_table_name}.tsv and reduce it to only the selected samples
    if [ -n "~{input_table}" ]; then
      cp "~{input_table}" "~{terra_table_name}.tsv"
    else
      bioforklift download "~{terra_table_name}" \
        --project "~{terra_project_name}" \
        --workspace "~{terra_workspace_name}" \
        --samples "${SAMPLE_IDS[@]}"
    fi

    # builds a metadata.csv (and FASTA) for every sample, validates each sample on its own with `seqsender prep`,
    # and combines the samples that pass into the 'submit' directory seqsender submits from.
    build_and_prep.py \
      --metadata_file "~{terra_table_name}.tsv" \
      --config_file "~{seqsender_config}" \
      --sample_ids "${SAMPLE_IDS[@]}" \
      ~{"--bioproject " + bioproject} \
      --assembly_fastas "${ASSEMBLY_FASTAS[@]}" \
      --read1_fastqs "${READ1_FASTQS[@]}" \
      --read2_fastqs "${READ2_FASTQS[@]}" \
      ~{"--sra_transfer_bucket " + sra_transfer_gcp_bucket} \
      ~{"--organism " + organism} \
      ~{if submit_to_biosample then '--biosample' else ''} \
      ~{if submit_to_sra then '--sra' else ''} \
      ~{if submit_to_genbank then '--genbank' else ''} \
      ~{if prep_only then '--prep_only' else ''}
  >>>
  output {
    File? combined_metadata = "submit/metadata.csv"
    File? combined_fasta = "submit/sequence.fasta"
    File? tarball = "seqsender_build_and_prep.tar.gz"
    File validation_summary = "validation_summary.tsv"
    Int num_samples_passed = read_int("NUM_PASSED")
    Int num_samples_failed = read_int("NUM_FAILED")
    File log = "seqsender_build_and_prep.log"
    String seqsender_version = read_string("VERSION")
    String seqsender_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 0
  }
}

task seqsender_submit {
  input {
    File submission_tarball # `tarball` (seqsender_build_and_prep.tar.gz) from seqsender_build_and_prep
    String seqsender_key # the key that encrypted the password in seqsender_config
    Boolean submit_to_production

    String organism
    Boolean submit_to_biosample
    Boolean submit_to_sra
    Boolean submit_to_genbank

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/seqsender:1.5.2-dev"
    Int cpu = 2
    Int memory = 4
    Int disk_size = 100
  }
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail

    # 'submit/' directory holds 'metadata.csv', 'sequence.fasta', 'config.yaml', and 'raw_reads/' when the reads are uploaded over FTP
    tar -xzf "~{submission_tarball}"

    NAME="seqsender_submit_$(date -u +%Y%m%d_%H%M%S)_~{organism}"
    echo "${NAME}" | tee SUBMISSION_NAME

    RETURN_CODE=0
    micromamba run -n seqsender seqsender.py submit \
      ~{if submit_to_biosample then '--biosample' else ''} \
      ~{if submit_to_sra then '--sra' else ''} \
      ~{if submit_to_genbank then '--genbank' else ''} \
      --organism "~{organism}" \
      --submission_dir submit \
      --submission_name "${NAME}" \
      --config_file submit/config.yaml \
      --metadata_file submit/metadata.csv \
      ~{if submit_to_genbank then '--fasta_file submit/sequence.fasta' else ''} \
      --key "~{seqsender_key}" \
      ~{if submit_to_production then '' else '--test'} \
      2>&1 | tee seqsender_submit.log || RETURN_CODE=$?

    # seqsender adds a row to submission_log.csv after each database is uploaded, so without it nothing was submitted
    if [ ! -f submit/submission_log.csv ]; then
      echo "ERROR: seqsender submit failed before anything was submitted; see seqsender_submit.log" >&2
      exit 1
    fi
    # a failure after some databases were uploaded still saves the state, so seqsender_poll can follow them
    if [ "${RETURN_CODE}" -eq 0 ]; then
      echo "SUBMITTED"
    else
      echo "ERROR"
    fi | tee SUBMIT_STATUS

    # the reads are only used by `seqsender submit`, so the `seqsender_state` for polling leaves them out
    tar -czf seqsender_state.tar.gz --exclude=submit/raw_reads submit
  >>>
  output {
    File tarball = "seqsender_state.tar.gz"
    String submission_name = read_string("SUBMISSION_NAME")
    String status = read_string("SUBMIT_STATUS")
    File log = "seqsender_submit.log"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 0
  }
}