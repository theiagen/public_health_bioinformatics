version 1.0

task bakta {
  input {
    File assembly
    String samplename
    String db_type = "light" # User choice for database type: "light" (default) or "full"
    Int cpu = 8
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bakta:1.10.3"
    File bakta_light_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_light_db_2024-01-20.tar.gz"
    File bakta_full_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/bakta_full_db_2024-01-20.tar.gz"
    Boolean proteins = false #  proteins: Fasta file of trusted protein sequences for CDS annotation
    Boolean compliant = false
    File? prodigal_tf # prodigal_tf: Prodigal training file to use for CDS prediction
    String? bakta_opts # bakta_opts: any additional bakta arguments
  }
  command <<<
  date | tee DATE
  bakta --version | tee BAKTA_VERSION
  
  # Determine which database to use
  if [[ "~{db_type}" == "light" ]]; then
    cp ~{bakta_light_db} db.tar.gz
  else
    cp ~{bakta_full_db} db.tar.gz
  fi

  # Extract Bakta DB
  mkdir db
  time tar xzvf db.tar.gz --strip-components=1 -C ./db

  # Install amrfinderplus db
  amrfinder_update --database db/amrfinderplus-db
  amrfinder --database_version | tee AMRFINDER_DATABASE_VERSION

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
  
  # rename gff3 to gff for compatibility with downstream analysis (pirate)
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
    maxRetries: 3
  }
}