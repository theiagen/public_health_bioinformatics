version 1.0

task submit_ena_data {
  input{
    # Input arguments
    String ena_username # ENA username
    String ena_password # ENA password

    File metadata_accessions # Metadata spreadsheet (TSV format) with ENA accessions

    String sample_id_column # Column name containing sample IDs
    Array[String] samples # Comma-separated list of sample IDs to include (default: all)

    # Column mappings if user wants to map Terra columns to ENA columns
    File? column_mappings # TSV file containing mappings from Terra columns to ENA columns

    # Column name configuration (used if column-mappings not provided)
    String read1_column = "read1" # Column name for read1/first fastq files (default: read1)
    String read2_column = "read2" # Column name for read2/second fastq files (default: read2)
    String bam_column = "bam_file" # Column name for BAM files (default: bam_file)
    String cram_column = "cram_file" # Column name for CRAM files (default: cram_file)

    # Pass in if not in the table, study_accession always required -- similar to how terra2_ncbi works
    String study_accession # ENA study accession
    String? experiment_name # Default experiment name template
    String library_strategy = "WGS" # Default library strategy (default: WGS)
    String library_source = "GENOMIC" # Default library source (default: GENOMIC)
    String library_selection = "RANDOM" # Default library selection (default: RANDOM)
    String platform = "ILLUMINA" # Default sequencing platform (default: ILLUMINA)
    String? instrument # Default sequencing instrument

    # Allow continuing even if some samples have missing metadata
    # Default is to fail if any samples are missing required metadata, unless this flag is set
    Boolean? allow_missing
    Boolean? test_submit

    Int disk_size = 100
    Int cpu = 1
    Int memory = 2
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra_to_ena:0.4"
  }
  command <<<
    set -euo pipefail

    # prepare metadata for ENA submission
    prepare_ena_data.py \
      ~{"--input " + metadata_accessions }  \
      --output "prepped_ena_data.tsv" \
      --file-paths "file_paths.json" \
      --excluded "excluded_samples.tsv" \
      ~{"--sample-id-column " + sample_id_column} \
      --samples ~{sep=',' samples} \
      ~{"--column-mappings " + column_mappings} \
      ~{"--read1-column " + read1_column} \
      ~{"--read2-column " + read2_column} \
      ~{"--bam-column " + bam_column} \
      ~{"--cram-column " + cram_column} \
      ~{"--study-accession " + study_accession} \
      ~{"--experiment-name " + experiment_name} \
      ~{"--library-strategy " + library_strategy} \
      ~{"--library-source " + library_source} \
      ~{"--library-selection " + library_selection} \
      ~{"--platform " + platform} \
      ~{"--instrument " + instrument} \
      ~{true="--allow-missing" false="" allow_missing}

    # localize files from a remote path (gs://) to a local path
    localize_files.py \
      --file-paths "file_paths.json"

    # Submit the data to ENA
    bulk_webincli.py \
      ~{"--username " + ena_username} \
      ~{"--password " + ena_password} \
      --webinCliPath /scripts/webin-cli.jar \
      --geneticContext "reads" \
      --spreadsheet "prepped_ena_data.tsv" \
      --directory . \
      --mode "submit" \
      ~{true="--test" false="" test_submit}

  >>>
  output {
    File prepped_ena_data = "prepped_ena_data.tsv"
    File file_paths_json = "file_paths.json"
    File excluded_samples = "excluded_samples.tsv"
    Array[File]? manifest_files = glob("manifests/Manifest_*.txt")
    Array[File]? manifest_log_err_files = glob("manifests/*-report/*.err")
    Array[File]? manifest_log_out_files = glob("manifests/*-report/*.out")
    File manifest_all_errors = "failed_validation.txt"
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