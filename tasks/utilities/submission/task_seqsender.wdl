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

    String seqsender_organism_tag
    Boolean submit_to_biosample
    Boolean submit_to_sra
    Boolean submit_to_genbank

    Boolean validation_only

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
      ~{"--organism " + seqsender_organism_tag} \
      ~{if submit_to_biosample then '--biosample' else ''} \
      ~{if submit_to_sra then '--sra' else ''} \
      ~{if submit_to_genbank then '--genbank' else ''} \
      ~{if validation_only then '--validation_only' else ''}
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

    String seqsender_organism_tag
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

    NAME="seqsender_submit_$(date -u +%Y%m%d_%H%M%S)_~{seqsender_organism_tag}"
    echo "${NAME}" | tee SUBMISSION_NAME

    RETURN_CODE=0
    micromamba run -n seqsender seqsender.py submit \
      ~{if submit_to_biosample then '--biosample' else ''} \
      ~{if submit_to_sra then '--sra' else ''} \
      ~{if submit_to_genbank then '--genbank' else ''} \
      --organism "~{seqsender_organism_tag}" \
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

    # the reads are only used by `seqsender submit`, so the checkpoint for polling leaves them out
    tar -czf seqsender_checkpoint.tar.gz --exclude=submit/raw_reads submit
  >>>
  output {
    File tarball = "seqsender_checkpoint.tar.gz"
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

task seqsender_poll_and_wait {
  input {
    File checkpoint_tarball # `tarball` from seqsender_submit or the newest seqsender_poll_and_wait
    String seqsender_key # the key that encrypted the password in seqsender_config
    Float max_wait_hours
    Float poll_interval_minutes

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/seqsender:1.5.2-dev"
    Int cpu = 1
    Int memory = 2
    Int disk_size = 10
  }
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail

    micromamba run -n seqsender \
    seqsender.py version | grep -o "[0-9][0-9.]*[0-9]" | head -n 1 | tee VERSION || true

    # 'submit/' directory holds the submission files and the submission_log.csv seqsender uses to track progress
    tar -xzf "~{checkpoint_tarball}"

    # checks NCBI ftp until every database is processed, one has an error, or max_wait_hours is reached
    RETURN_CODE=0
    poll_and_wait.py \
      --submit_dir submit \
      --key "~{seqsender_key}" \
      --max_wait_hours ~{max_wait_hours} \
      --poll_interval_minutes ~{poll_interval_minutes} || RETURN_CODE=$?

    # saved even if polling stopped early
    tar -czf seqsender_checkpoint.tar.gz submit
    exit "${RETURN_CODE}"
  >>>
  output {
    File tarball = "seqsender_checkpoint.tar.gz"
    String submission_name = read_string("SUBMISSION_NAME")
    String final_status = read_string("FINAL_STATUS")
    String database_statuses = read_string("DATABASE_STATUSES")
    File log = "seqsender_poll_and_wait.log"
    File submission_metadata = "submit/metadata.csv"
    File? submission_fasta = "submit/sequence.fasta"
    File submission_log_csv = "submit/submission_log.csv"
    File? submission_status_report = "submit/~{submission_name}/submission_files/submission_status_report.csv"
    File? report_xml_biosample = "submit/~{submission_name}/submission_files/BIOSAMPLE/report.xml"
    File? report_xml_sra = "submit/~{submission_name}/submission_files/SRA/report.xml"
    File? report_xml_genbank = "submit/~{submission_name}/submission_files/GENBANK/report.xml"
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

task seqsender_upload_accessions {
  input {
    String terra_project_name
    String terra_workspace_name

    File submission_status_report # `submission_status_report` from seqsender_poll_and_wait
    File submission_metadata # `submission_metadata` from seqsender_poll_and_wait; links each sample name to its Terra row

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/seqsender:1.5.2-dev"
    Int cpu = 1
    Int memory = 2
    Int disk_size = 10
  }
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail

    # adds each database's accession column on the submitted rows only; all other columns/rows untouched
    python3 <<CODE
    import json
    import pandas as pd
    from bioforklift.terra import Terra

    accession_columns = {
      "bs-sample_name": "biosample_accession",
      "sra-sample_name": "sra_accession",
      "gb-sample_name": "genbank_accession"
    }

    status_df = pd.read_csv("~{submission_status_report}", dtype=str, keep_default_na=False)
    metadata_df = pd.read_csv("~{submission_metadata}", dtype=str, keep_default_na=False)

    # the 'metadata.csv' file used for submission keeps the original Terra table's 'entity:{table_name}_id' column
    entity_column = next(column for column in metadata_df.columns if column.startswith("entity:"))
    table_name = entity_column.removeprefix("entity:").removesuffix("_id")

    terra = Terra(
      source_workspace="~{terra_workspace_name}",
      source_project="~{terra_project_name}",
      destination_workspace="~{terra_workspace_name}",
      destination_project="~{terra_project_name}"
    )

    # one row per submitted sample; database accession headers will be added where applicable
    upload_df = metadata_df[[entity_column]].drop_duplicates()

    # seqsender will only list accession headers for databases it attempts to submit to
    for name_column, accession_column in accession_columns.items():
      if name_column not in status_df:
        continue

      # flu samples have one row per segment; the BioSample and SRA accessions will repeat; GenBank accessions should be unique
      df = status_df.loc[status_df[accession_column] != "", [name_column, accession_column]].drop_duplicates()
      df = df.merge(metadata_df[[entity_column, name_column]].drop_duplicates(), on=name_column)
      if df.empty:
        print(f"No {accession_column} values to upload")
        continue

      if accession_column == "genbank_accession":
        # one JSON map of "gb-sample_name": "genbank_accession" per sample
        accessions_df = (
          df.groupby(entity_column, sort=False)[[name_column, accession_column]]
          .apply(lambda records: json.dumps(dict(records.values)))
        )
      else:
        # one BioSample/SRA accession per sample
        accessions_df = df.set_index(entity_column)[accession_column]
      upload_df[accession_column] = upload_df[entity_column].map(accessions_df)
      print(f"Found {len(accessions_df)} {accession_column} values")

    # a sample missing an accession gets a blank cell, which Terra skips rather than clearing the existing value
    if len(upload_df.columns) == 1:
      print("No accessions to upload")
    else:
      terra.entities.upload_entities(data=upload_df, target=table_name)
      print(f"Uploaded {', '.join(upload_df.columns[1:])} for {len(upload_df)} samples to {table_name}")
    CODE
  >>>
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