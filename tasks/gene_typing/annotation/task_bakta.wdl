version 1.0

task bakta {
  input {
    File assembly
    String samplename
    Int cpu = 8
    Int memory = 24 
    Int disk_size = 200
    File bakta_db_selected
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bakta:1.10.3"
    Boolean compliant = false
    File? prodigal_tf # Prodigal_tf: Prodigal training file to use for CDS prediction
    File? proteins # Proteins: Fasta file of trusted protein sequences for CDS annotation
    String? bakta_opts # Additional Bakta arguments
  }
  command <<<  
  set -euo pipefail

  date | tee DATE
  bakta --version | tee BAKTA_VERSION

  # Extract Bakta DB
  mkdir -p db
  tar xzvf ~{bakta_db_selected} --strip-components=1 -C ./db

  # Run Bakta
  bakta \
    ~{bakta_opts} \
    --db db/ \
    --threads ~{cpu} \
    --prefix ~{samplename} \
    --output ~{samplename} \
    ~{true='--compliant' false='' compliant} \
    ~{'--proteins ' + proteins} \
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
    File bakta_plot = "~{samplename}/~{samplename}.png"
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
