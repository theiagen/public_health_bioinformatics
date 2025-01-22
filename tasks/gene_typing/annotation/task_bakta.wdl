version 1.0

task bakta {
  input {
    File assembly
    String samplename
    String db_type = "light" # User choice for database type: light (default) or full
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
    String docker = "bakta:1.10.3-local"
    String bakta_light_db_url = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_light_db_2024-01-20.tar.gz"
    String bakta_full_db_url = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_full_db_2024-01-20.tar.gz"
    Boolean proteins = false # Proteins: Fasta file of trusted protein sequences for CDS annotation
    Boolean compliant = false
    File? prodigal_tf # Prodigal_tf: Prodigal training file to use for CDS prediction
    String? bakta_opts # Additional Bakta arguments
  }
  command <<<  
  set -euo pipefail
  date | tee DATE
  bakta --version | tee BAKTA_VERSION

   # Define database URLs
  BAKTA_LIGHT_DB_URL="gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_light_db_2024-01-20.tar.gz"
  BAKTA_FULL_DB_URL="gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_full_db_2024-01-20.tar.gz"

  # Debug statement for database type
  echo "Using database type: ~{db_type}" | tee DB_TYPE

  # Determine database URL
  if [[ "~{db_type}" == "light" ]]; then
    echo "Using light database" | tee DB_DOWNLOAD_LOG
    DB_URL="$BAKTA_LIGHT_DB_URL"
  else
    echo "Using full database" | tee DB_DOWNLOAD_LOG
    DB_URL="$BAKTA_FULL_DB_URL"
  fi

  # Download the selected database
  echo "Downloading database from: $DB_URL" | tee -a DB_DOWNLOAD_LOG
  gcloud storage cp $DB_URL db.tar.gz


  # Extract Bakta DB
  mkdir db
  time tar xzvf db.tar.gz --strip-components=1 -C ./db

  # Run Bakta
  bakta \
    ~{bakta_opts} \
    --db db/ \
    --threads ~{cpu} \
    --prefix ~{samplename} \
    --output ~{samplename} \
    ~{true='--compliant' false='' compliant} \
    ~{true='--proteins' false='' proteins} \
    ~{'--prodigal-tf ' + prodigal_tf} \
    ~{assembly}

  # Rename gff3 to gff for compatibility with downstream analysis
  mv "~{samplename}/~{samplename}.gff3" "~{samplename}/~{samplename}.gff"
  >>>
  output {
    File bakta_embl = "~{samplename}/~{samplename}.embl"
    File bakta_faa = "~{samplename}/~{samplename}.faa"
    File bakta_ffn = "~{samplename}/~{samplename}.ffn"
    File bakta_fna = "~{samplename}/~{samplename}.fna"
    File bakta_gbff = "~{samplename}/~{samplename}.gbff"
    File bakta_gff3 = "~{samplename}/~{samplename}.gff"
    File bakta_hypotheticals_faa = "~{samplename}/~{samplename}.hypotheticals.faa"
    File bakta_hypotheticals_tsv = "~{samplename}/~{samplename}.hypotheticals.tsv"
    File bakta_tsv = "~{samplename}/~{samplename}.tsv"
    File bakta_txt = "~{samplename}/~{samplename}.txt"
    String bakta_version = read_string("BAKTA_VERSION")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
  }
}
