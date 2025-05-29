version 1.0

task register_ena_samples {
  input {
    File metadata # Metadata spreadsheet (TSV format)
    String ena_username
    String ena_password
    String study_accession # ENA study accession
    String sample_id_column # Column name containing sample IDs
    String sample_type # must be "prokaryotic_pathogen" or "virus_pathogen"

    String? center # Center name for submission
    File? column_mappings # TSV file with column name mappings from Terra to ENA
    Boolean? allow_missing # Allow missing metadata for some samples

    Int batch_size = 100 # Number of samples to submit in one batch (default: 100)
    Int disk_size = 100
    Int cpu = 1
    Int memory = 2
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra_to_ena:0.5"
  }
  command <<<
      set -euo pipefail

      register_samples_ena.py \
        ~{"--metadata " + metadata} \
        ~{"--username " + ena_username} \
        ~{"--password " + ena_password} \
        ~{"--study " + study_accession} \
        ~{"--sample-id-column " + sample_id_column} \
        ~{"--sample-type " + sample_type} \
        ~{"--center " + center} \
        ~{true="--allow-missing" false="" allow_missing} \
        ~{"--column-mappings " + column_mappings} \
        ~{"--batch-size " + batch_size} \
        --output .

  >>>
  output {
      File accessions = "accessions.tsv"
      File metadata_accessions = "metadata_with_accessions.tsv"
      File registration_summary = "submission_summary.txt"
      File registration_log = "submission.log"
      String registration_success = read_string("success.txt")
      String docker_image = docker
  }
  runtime {
    cpu: cpu
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    docker: docker
  }
}