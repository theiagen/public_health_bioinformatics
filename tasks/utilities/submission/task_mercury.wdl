version 1.0

task mercury {
  input {
    File data_table
    String table_name
    Array[String] samplenames
    String gcp_bucket_uri
    String terra_project_name
    String terra_workspace_name

    # optional parameters
    String organism = "sars-cov-2"
    String output_name = "mercury"
    Boolean skip_county = false
    Boolean skip_ncbi = false
    Boolean using_clearlabs_data = false
    Boolean using_reads_dehosted = false
    Boolean usa_territory = false
    Boolean single_end = false
    Int vadr_alert_limit = 0
    Int number_N_threshold = 5000
    String authors = ""
    String bioproject_accession = ""
    String continent = ""
    String country = ""
    String host_disease = ""
    String isolation_source = ""
    String library_selection = ""
    String library_source = ""
    String library_strategy = ""
    String purpose_of_sequencing = ""
    String state = ""
    String submitting_lab = ""
    String submitting_lab_address = ""
    String amplicon_primer_scheme = ""
    String amplicon_size = ""
    String instrument_model = ""
    String library_layout = ""
    String seq_platform = ""
    String gisaid_submitter = ""
    String submitter_email = ""
    String metadata_organism = ""
  
    # runtime parameters
    Int cpu = 2
    Int disk_size = 100
    Int memory = 8
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/mercury:1.1.3"
  }
  meta {
    volatile: true
  }
  command <<<
    #set -euo pipefail to avoid silent failure
    set -euo pipefail

    python3 /mercury/mercury/mercury.py -v | tee VERSION

    python3 /mercury/mercury/mercury.py \
      ~{data_table} "~{table_name}_id" ~{sep=',' samplenames} \
      --gcp_bucket_uri ~{gcp_bucket_uri} \
      ~{"--organism " + organism} \
      ~{"--output_prefix " + output_name} \
      ~{true="--skip_county" false="" skip_county} \
      ~{true="--skip_ncbi" false="" skip_ncbi} \
      ~{true="--using_clearlabs_data" false="" using_clearlabs_data} \
      ~{true="--using_reads_dehosted" false="" using_reads_dehosted} \
      ~{true="--usa_territory" false="" usa_territory} \
      ~{true="--single_end" false="" single_end} \
      ~{"--vadr_alert_limit " + vadr_alert_limit} \
      ~{"--number_n_threshold " + number_N_threshold} \
      ~{"--authors '" + authors + "'"} \
      ~{"--bioproject_accession '" + bioproject_accession + "'"} \
      ~{"--continent '" + continent + "'"} \
      ~{"--country '" + country + "'"} \
      ~{"--host_disease '" + host_disease + "'"} \
      ~{"--isolation_source '" + isolation_source + "'"} \
      ~{"--library_selection '" + library_selection + "'"} \
      ~{"--library_source '" + library_source + "'"} \
      ~{"--library_strategy '" + library_strategy + "'"} \
      ~{"--purpose_of_sequencing '" + purpose_of_sequencing + "'"} \
      ~{"--state '" + state + "'"} \
      ~{"--submitting_lab '" + submitting_lab + "'"} \
      ~{"--submitting_lab_address '" + submitting_lab_address + "'"} \
      ~{"--amplicon_primer_scheme '" + amplicon_primer_scheme + "'"} \
      ~{"--amplicon_size '" + amplicon_size + "'"} \
      ~{"--instrument_model '" + instrument_model + "'"} \
      ~{"--library_layout '" + library_layout + "'"} \
      ~{"--seq_platform '" + seq_platform + "'"} \
      ~{"--gisaid_submitter '" + gisaid_submitter + "'"} \
      ~{"--submitter_email '" + submitter_email + "'"} \
      ~{"--metadata_organism '" + metadata_organism + "'"} \
      --debug

    # write out excluded samples file to the stdout
    cat ~{output_name}_excluded_samples.tsv

    # output to the initial Terra table with the updated metadata
    python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{terra_project_name}" --workspace "~{terra_workspace_name}" --tsv ~{output_name}_terra_table_to_upload.tsv

    # provide a lowercased organism variable for use later
    echo "~{organism}" | tr '[:upper:]' '[:lower:]' > ORGANISM_NAME

    # change bankit_metadata to tsv if exists
    if [ -f ~{output_name}.src ]; then
      mv ~{output_name}.src ~{output_name}.tsv
    fi
  >>>
  output {
    String mercury_version = read_string("VERSION")
    String organism_name = read_string("ORGANISM_NAME")
    File? bankit_fasta = "~{output_name}_bankit_combined.fasta"
    File? bankit_metadata = "~{output_name}.tsv"
    File? biosample_metadata = "~{output_name}_biosample_metadata.tsv"
    File? excluded_samples = "~{output_name}_excluded_samples.tsv"
    File? genbank_fasta = "~{output_name}_genbank_untrimmed_combined.fasta"
    File? genbank_metadata = "~{output_name}_genbank_metadata.tsv"
    File? gisaid_fasta = "~{output_name}_gisaid_combined.fasta"
    File? gisaid_metadata = "~{output_name}_gisaid_metadata.csv"
    File? sra_metadata = "~{output_name}_sra_metadata.tsv"
    File? terra_table = "~{output_name}_terra_table_to_upload.tsv"
  }
  runtime {
    cpu: cpu
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: docker
    memory: memory
    preemptible: 0 # this task has may have a long time and shouldn't be preempted
  }
}