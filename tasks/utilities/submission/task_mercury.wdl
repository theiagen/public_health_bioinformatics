version 1.0

task mercury {
  input {
    File data_table
    String table_name
    Array[String] samplenames
    String gcp_bucket_uri

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
    
    # runtime parameters
    Int cpu = 2
    Int disk_size = 100
    Int memory = 8
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/mercury:1.0.7"
  }
  meta {
    volatile: true
  }
  command <<<
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
      --debug

    # write out excluded samples file to the stdout
    cat ~{output_name}_excluded_samples.tsv
  >>>
  output {
    String mercury_version = read_string("VERSION")
    File? bankit_fasta = "~{output_name}_bankit_combined.fasta"
    File? bankit_metadata = "~{output_name}.src"
    File? biosample_metadata = "~{output_name}_biosample_metadata.tsv"
    File? excluded_samples = "~{output_name}_excluded_samples.tsv"
    File? genbank_fasta = "~{output_name}_genbank_untrimmed_combined.fasta"
    File? genbank_metadata = "~{output_name}_genbank_metadata.tsv"
    File? gisaid_fasta = "~{output_name}_gisaid_combined.fasta"
    File? gisaid_metadata = "~{output_name}_gisaid_metadata.csv"
    File? sra_metadata = "~{output_name}_sra_metadata.tsv"
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