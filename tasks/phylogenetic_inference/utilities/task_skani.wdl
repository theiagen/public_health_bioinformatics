version 1.0

task skani {
  input{
    File assembly_fasta
    String samplename
    #File skani_db = "gs://theiagen-large-public-files-rp/terra/databases/skani/skani_20250314.tar.gz"
    File skani_db
    Int disk_size = 100
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/skani:0.2.2"
  }
  command <<<
    set -euo pipefail

    # extract skani db
    tar -xf ~{skani_db}
    untarred_skani_db=$(basename ~{skani_db} .tar)

    # get version
    skani --version | tee VERSION

    skani search \
      -d ${untarred_skani_db} \
      -q ~{assembly_fasta} \
      -s 50 \
      --no-learned-ani \
      --robust \
      --detailed \
      --ci \
      --no-marker-index \
      -o ~{samplename}_skani_results.tsv

    # add new column header
    new_header=$(awk 'NR==1 {OFS="\t"; print $0, "ANI_x_Total_bases"}' ~{samplename}_skani_results.tsv)

    # create a new column and sort by the product of ANI (col 3) and Total bases covered (col 20)
    awk -F'\t' 'NR > 1 {OFS="\t"; new_col = sprintf("%.2f", $3 * $20); print $0, new_col}' ~{samplename}_skani_results.tsv | \
      sort -t$'\t' -k21,21nr | \
      { echo "$new_header"; cat -; } > ~{samplename}_skani_results_sorted.tsv

    # get the file name for the top hit in sorted skani results: 1st column, 2nd line in file
    top_hit_name=$(basename $(awk 'NR == 2 {print $1}' ~{samplename}_skani_results_sorted.tsv))

    if [ -z "$top_hit_name" ]; then
      echo "ERROR: No hits found in skani results"
      exit 1
    else
      # get accession number from file name (and version number if it exists)
      echo $top_hit_name | awk -F'[.]' '{print $1 ($2 ~ /^[0-9]+$/ ? "."$2 : "")}' | tee TOP_ANI_ACCESSION
    fi

  >>>
  output{
    File skani_report = "~{samplename}_skani_results_sorted.tsv"
    String skani_top_ani_accession = read_string("TOP_ANI_ACCESSION")
    String skani_database = skani_db
    String skani_version = read_string("VERSION")
    String skani_docker = docker
  }
  runtime{
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}