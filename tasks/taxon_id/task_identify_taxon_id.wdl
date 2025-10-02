version 1.0

task identify_taxon_id {
  input {
    String taxon # can be taxon id (int) or organism name (string)
    String? rank # limit input taxon to the user specified rank
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

    echo "DEBUG: Obtaining taxon report for taxon: ~{taxon} and rank: ~{rank}"
    # Note "jq -c '.reports[]'" acts identically to --as-json-lines
    datasets summary taxonomy taxon ~{'"' + taxon + '"'} | \
    tee ncbi_taxon_summary.json | \
    jq -c '.reports[]' | \
    dataformat tsv taxonomy --template tax-summary > ncbi_taxon_summary.tsv

    python3 <<CODE
    import json
    import sys

    with open("ncbi_taxon_summary.json") as f:
      data = json.load(f)
      report_list = data['reports']

      for report in report_list:
        query_taxon = ', '. join([q for q in report['query']])
        query_rank = "~{rank}".lower()
        report = report['taxonomy']

        # Raw taxon name, id, and rank is based on the original user input query
        raw_taxon_name = report['current_scientific_name'].get('name', '')
        raw_taxon_id = report.get('tax_id', '')
        if 'rank' in report:
          raw_taxon_rank = report['rank'].lower()
        # Taxon is below species level. Set to 'no rank' (complicit with NCBI Taxonomy database conventions).
        elif 'species' in report['classification']:
          raw_taxon_rank = 'no rank'
        else:
          raw_taxon_rank = 'N/A'

        # Reported taxon name, id, and rank is based on the ranked user input query (if provided/found)
        # if no rank provided, default to raw taxon rank, unless taxon is below species level (no rank) then set to "species"
        reported_taxon_rank = query_rank if query_rank else ('species' if raw_taxon_rank == 'no rank' else raw_taxon_rank)
        if query_rank and (reported_taxon_rank not in report['classification']):
          reported_taxon_rank = reported_taxon_name = reported_taxon_id = 'N/A'
          print(f"ERROR: Input taxon rank '{query_rank}' is not valid (too specific) for taxon: '{query_taxon}'.")
          sys.exit(1)
        else:
          reported_taxon_name = report['classification'].get(reported_taxon_rank, {}).get('name', '')
          reported_taxon_id = report['classification'].get(reported_taxon_rank, {}).get('id', '')

    outputs = {
      "TAXON_ID": str(reported_taxon_id),
      "TAXON_NAME": str(reported_taxon_name),
      "TAXON_RANK": str(reported_taxon_rank),
      "RAW_TAXON_ID": str(raw_taxon_id),
      "RAW_TAXON_NAME": str(raw_taxon_name),
      "RAW_TAXON_RANK": str(raw_taxon_rank)
    }
    for filename, value in outputs.items():
      with open(filename, "w") as f:
        f.write(value)
    CODE
    echo "DEBUG: Reported TAXON_ID: $(cat TAXON_ID), TAXON_NAME: $(cat TAXON_NAME), TAXON_RANK: $(cat TAXON_RANK)"

    if [ ! -s "ncbi_taxon_summary.tsv" ]; then
      echo "ERROR: no taxon summary found for taxon: ~{taxon} and rank: ~{rank}"
      exit 1
    fi

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
      }' ncbi_genome_summary.tsv | grep -Po "\d+" > AVG_GENOME_LENGTH
    fi
  >>>
  output {
    File taxon_summary_tsv = "ncbi_taxon_summary.tsv"
    File taxon_summary_json = "ncbi_taxon_summary.json"
    File genome_summary_tsv = "ncbi_genome_summary.tsv"
    String taxon_id = read_string("TAXON_ID")
    String taxon_name = read_string("TAXON_NAME")
    String taxon_rank = read_string("TAXON_RANK")
    String raw_taxon_id = read_string("RAW_TAXON_ID")
    String raw_taxon_name = read_string("RAW_TAXON_NAME")
    String raw_taxon_rank = read_string("RAW_TAXON_RANK")
    String avg_genome_length = read_string("AVG_GENOME_LENGTH")
    String ncbi_datasets_accession = read_string("NCBI_ACCESSION")
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