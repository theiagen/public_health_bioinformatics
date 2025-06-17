version 1.0

task identify_taxon_id {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    String? rank # limit input taxon to the user specified rank
    Int summary_limit = 100 # limit the number of genomes to summarize
    Boolean use_ncbi_virus = false
    Boolean complete = true
    Boolean refseq = true
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
    datasets summary taxonomy taxon ~{'"' + taxon + '"'} \
      ~{"--rank " + rank} \
      --as-json-lines | \
    dataformat tsv taxonomy \
      --template tax-summary > ncbi_taxon_summary.tsv

    # check if the taxon summary file is empty
    if [ ! -s "ncbi_taxon_summary.tsv" ]; then
      echo "ERROR: no taxon summary found for taxon: ~{taxon} and rank: ~{rank}"
      exit 1
    else
      # check if ncbi_taxon_summary.tsv is longer than 2 lines (header + 1 data line)
      # if so, then the taxon rank (either calculated or provided) is too specific and we need to exit with an error
      if [ $(wc -l < ncbi_taxon_summary.tsv) -gt 2 ]; then
        echo "ERROR: input taxon rank '~{rank}' is not valid (too specific) for taxon: '~{taxon}'." 1>&2
        exit 1
      fi
      # extract the taxid ($2) tax name ($3) and rank ($5) from the output tsv file
      # skip the header line. if the input taxon rank is invalid (too specific), there will be more than one line listed here.
      awk -F'\t' 'NR == 2 {print $2}' ncbi_taxon_summary.tsv > TAXON_ID
      awk -F'\t' 'NR == 2 {print $3}' ncbi_taxon_summary.tsv > TAXON_NAME
      awk -F'\t' 'NR == 2 {print $5}' ncbi_taxon_summary.tsv > TAXON_RANK
    fi

    # get list of {summary_limit} genomes from the specified taxon and calculate average genome length
    if ~{use_ncbi_virus}; then
      # virus-genome --fields: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/using-dataformat/virus-data-reports/
      datasets summary \
        virus \
        genome taxon ~{'"' + taxon + '"'} \
        ~{"--limit " + summary_limit} \
        ~{true="--complete-only" false="" complete} \
        ~{true="--refseq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv virus-genome \
        --fields accession,completeness,is-annotated,length,sourcedb,virus-name,virus-tax-id \
         > ncbi_genome_summary.tsv
    else
      # genome --fields: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/using-dataformat/genome-data-reports/
      datasets summary \
        genome taxon ~{'"' + taxon + '"'} \
        ~{"--limit " + summary_limit} \
        ~{true="--assembly-level complete" false="" complete} \
        ~{true="--assembly-source RefSeq" false="" refseq} \
        --as-json-lines | \
      dataformat tsv genome \
        --fields accession,assmstats-number-of-contigs,checkm-completeness,assmstats-total-sequence-len	\
         > ncbi_genome_summary.tsv
    fi

    # report if the file is empty
    if [ ! -s "ncbi_genome_summary.tsv" ]; then
      echo "ERROR: no assemblies found for taxon: ~{taxon}"
      echo "N/A" > NCBI_ACCESSION
      echo "0" > AVG_GENOME_LENGTH
    else
      ncbi_accession=$(awk -F'\t' 'NR == 2 {print $1}' ncbi_genome_summary.tsv)
      echo "DEBUG: ncbi_accession: ${ncbi_accession}"
      echo $ncbi_accession > NCBI_ACCESSION

      # take the average of of all genome lengths across each row
      awk -F'\t' '{
        if (NR > 1) { sum += $4; count++ }
      }
      END {
        if (count) { printf "%.0f", sum / count } else { print "No matching taxon id found" }
      }' ncbi_genome_summary.tsv > AVG_GENOME_LENGTH
    fi
  >>>
  output {
    File taxon_summary_tsv = "ncbi_taxon_summary.tsv"
    File genome_summary_tsv = "ncbi_genome_summary.tsv"
    String? taxon_id = read_string("TAXON_ID")
    String? taxon_name = read_string("TAXON_NAME")
    String? taxon_rank = read_string("TAXON_RANK")
    Int avg_genome_length = read_int("AVG_GENOME_LENGTH")
    String? ncbi_datasets_accession = read_string("NCBI_ACCESSION")
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