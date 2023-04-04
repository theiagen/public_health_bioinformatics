version 1.0

task ncbi_datasets_download_genome_accession {
  input {
    String ncbi_accession
    Int cpu = 1
    Int memory = 4
    String docker = "staphb/ncbi-datasets:14.13.2" # not the latest version, but it's hard to keep up w/ the frequent releases
    Int disk_size = 50
  }
  command <<<
  date | tee DATE
  datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

  # download FASTA file using ncbi_accession
  # '--assembly-version latest' ensures the most recent version is downloaded, not previous versions
  # NOTE: I have noticed that occasionally the same command may fail one moment with the error below and succeed the second time. usually.
  ### "Error: No assemblies found that match selection"
  datasets download genome accession \
    ~{ncbi_accession} \
    --filename ~{ncbi_accession}.zip \
    --assembly-version latest

  # unzip the archive and copy FASTA and JSON to PWD, rename in the process so output filenames are predictable 
  unzip ~{ncbi_accession}.zip
  cp -v ncbi_dataset/data/~{ncbi_accession}*/~{ncbi_accession}*.fna ./~{ncbi_accession}.fasta
  cp -v ncbi_dataset/data/assembly_data_report.jsonl ./~{ncbi_accession}.data_report.jsonl

  >>>
  output {
    File ncbi_datasets_assembly_fasta = "~{ncbi_accession}.fasta"
    File ncbi_datasets_assembly_data_report_json = "~{ncbi_accession}.data_report.jsonl"
    String ncbi_datasets_version = read_string("DATASETS_VERSION")
    String ncbi_datasets_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}