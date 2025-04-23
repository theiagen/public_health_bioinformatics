version 1.0

task identify_taxon_id {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    String? rank # limit input taxon to the user specified rank
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:16.38.1"
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    # fail hard
    set -euo pipefail

    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    echo "DEBUG: Obtaining taxon report for taxon: ~{taxon} and rank: ~{rank}"
    datasets summary taxonomy taxon ~{'"' + taxon + '"'} \ # <- im sorry
      ~{"--rank " + rank} \
      --as-json-lines | \
    dataformat tsv taxonomy \
      --template tax-summary > ncbi_taxon_summary.tsv

    # extract the taxid ($2) tax name ($3) and rank ($5) from the output tsv file
    # skip the header line. if the input taxon rank is invalid (too specific), there will be more than one line listed here.
    awk -F'\t' 'NR == 2 {print $2}' ncbi_taxon_summary.tsv > TAXON_ID
    awk -F'\t' 'NR == 2 {print $3}' ncbi_taxon_summary.tsv > TAXON_NAME
    awk -F'\t' 'NR == 2 {print $5}' ncbi_taxon_summary.tsv > TAXON_RANK

    # check if ncbi_taxon_summary.tsv is longer than 2 lines (header + 1 data line)
    # if so, then the taxon rank (either calculated or provided) is too specific and we need to exit with an error
    if [ $(wc -l < ncbi_taxon_summary.tsv) -gt 2 ]; then
      echo "ERROR: input taxon rank '~{rank}' is not valid (too specific) for taxon: '~{taxon}'." 1>&2
      exit 1
    fi
  >>>
  output {
    String taxon_id = read_string("TAXON_ID")
    String taxon_name = read_string("TAXON_NAME")
    String taxon_rank = read_string("TAXON_RANK")
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