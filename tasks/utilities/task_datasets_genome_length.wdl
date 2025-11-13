version 1.0

task datasets_genome_length {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    Int summary_limit = 100 # limit the number of genomes to summarize
    Boolean use_ncbi_virus = false
    Boolean complete = true
    Boolean refseq = true
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ncbi-datasets:18.9.0-python-jq"
    Int cpu = 1
    Int memory = 4
    Int disk_size = 50
  }
  command <<<
    # fail hard
    set -euo pipefail

    date | tee DATE
    datasets --version | sed 's|datasets version: ||' | tee DATASETS_VERSION

    echo "DEBUG: Generating genome summary for taxon: ~{taxon}"
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
      echo "0" > AVG_GENOME_LENGTH
    else
      # take the average of of all genome lengths across each row
      awk -F'\t' '{
        if (NR > 1) { sum += $4; count++ }
      }
      END {
        if (count) { printf "%.0f", sum / count } else { print "No matching taxon id found" }
      }' ncbi_genome_summary.tsv | grep -Po "\d+" > AVG_GENOME_LENGTH
    fi
  >>>
  output {
    File genome_summary_tsv = "ncbi_genome_summary.tsv"
    String avg_genome_length = read_string("AVG_GENOME_LENGTH")
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