version 1.0

task skani {
  input{
    File assembly_fasta
    String samplename
    File skani_db = "gs://theiagen-public-resources-rp/reference_data/databases/skani/skani_db_20251103.tar"
    String fasta_dir = "gs://theiagen-public-resources-rp/reference_data/databases/skani/viral_fna_20251103/"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/skani:0.2.2"
  }
  command <<<
    set -euo pipefail

    # extract skani db
    echo "DEBUG: Decompressing Skani DB"
    tar -xf ~{skani_db}
    untarred_skani_db=$(basename ~{skani_db} .tar)

    # get version
    skani --version | tee VERSION

    # set to PASS by default
    echo "PASS" | tee SKANI_STATUS

    # check if the fasta has contigs greater than 500 bp (otherwise skani will fail)
    echo "DEBUG: Checking assembly for contigs > 500 bp"
    sed -E 's/^>([^ ]+).*$/>\1/' ~{assembly_fasta} | \
      awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID"_extracted.fasta"} {print >> F}'
    compatible_fasta=0
    for fasta in *_extracted.fasta; do
      fa_length=$(grep -v '^>' $fasta | tr -d '\n' | wc -m)
      if [[ $fa_length -gt 500 ]]; then
        compatible_fasta=1
        echo "DEBUG: fasta has contigs > 500 bp"
        touch SKANI_WARNING
        break
      fi
    done

    # concatenate all contigs to spoof skani
    if [[ $compatible_fasta -eq 0 ]]; then
      echo "DEBUG: No contigs greater than 500 bp; concatenating supercontig for search"
      echo "All contigs < 500 bp" > SKANI_WARNING
      echo ">concatenated_supercontig" > concatenated.fasta
      for fasta in *_extracted.fasta; do
        grep -v '^>' $fasta >> concatenated.fasta
      done
      assembly_fasta="concatenated.fasta"

      # check if supercontig is greater than 500 bp
      fa_length=$(grep -v '^>' concatenated.fasta | tr -d '\n' | wc -m)
      if [[ $fa_length -gt 500 ]]; then
        echo "DEBUG: The combined length of the contigs exceeds 500 bp. (${fa_length} bp)"
      else
        echo "ERROR: The combined length of the contigs is less than 500 bp. (${fa_length} bp)"
        echo "FAIL" | tee SKANI_STATUS
      fi
    else
      assembly_fasta=~{assembly_fasta}
    fi

    # run skani
    echo "DEBUG: Running Skani"
    skani search \
      -d ${untarred_skani_db} \
      -q ${assembly_fasta} \
      -s 50 \
      --no-learned-ani \
      --robust \
      --detailed \
      --ci \
      --no-marker-index \
      -o ~{samplename}_skani_results.tsv

    # add new column header
    echo "DEBUG: Extracting Skani results"
    new_header=$(awk 'NR==1 {OFS="\t"; print $0, "ANI_x_Query_Coverage"}' ~{samplename}_skani_results.tsv)

    # initialize output files
    echo "N/A" > TOP_ACCESSION
    echo 0 > TOP_ANI
    echo 0 > TOP_QUERY_COVERAGE
    echo 0 > TOP_SCORE

    # create a new column and sort by the product of ANI (col 3) and Align_fraction_query (col 5)
    awk -F'\t' 'NR > 1 {OFS="\t"; new_col = sprintf("%.5f", $3 * $5 / 100); print $0, new_col}' ~{samplename}_skani_results.tsv | \
      sort -t$'\t' -k21,21nr | \
      { echo "$new_header"; cat -; } > ~{samplename}_skani_results_sorted.tsv

    if [[ $(wc -l < ~{samplename}_skani_results_sorted.tsv) -lt 2 ]]; then
      echo "ERROR: No hits found in skani results"
      echo "FAIL" | tee SKANI_STATUS
    else
      # get the file name for the top hit in sorted skani results: 1st column, 2nd line in file
      top_hit_name=$(basename $(awk 'NR == 2 {print $1}' ~{samplename}_skani_results_sorted.tsv))

      # get accession number from file name (and version number if it exists)
      echo $top_hit_name | awk -F'[.]' '{print $1 ($2 ~ /^[0-9]+$/ ? "."$2 : "")}' | tee TOP_ACCESSION
      head -n 2 ~{samplename}_skani_results_sorted.tsv | tail -n 1 | cut -f 3 | tee TOP_ANI
      head -n 2 ~{samplename}_skani_results_sorted.tsv | tail -n 1 | cut -f 5 | tee TOP_QUERY_COVERAGE
      head -n 2 ~{samplename}_skani_results_sorted.tsv | tail -n 1 | cut -f 21 | tee TOP_SCORE
    fi

  cat ~{fasta_dir}/ $(cat TOP_ACCESSION).fna > TOP_ASSEMBLY
  >>>
  output{
    File skani_report = "~{samplename}_skani_results_sorted.tsv"
    String skani_top_accession = read_string("TOP_ACCESSION")
    Float skani_top_ani = read_float("TOP_ANI")
    Float skani_top_query_coverage = read_float("TOP_QUERY_COVERAGE")
    Float skani_top_score = read_float("TOP_SCORE")
    String skani_reference_assembly = read_string("TOP_ASSEMBLY")
    String skani_database = skani_db
    String skani_warning = read_string("SKANI_WARNING")
    String skani_status = read_string("SKANI_STATUS")
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
