version 1.0

task checkv {
  meta {
    description: "Run CheckV on viral assemblies"
  }
  input {
    File assembly
    String samplename
    File checkv_db = "gs://theiagen-large-public-files-rp/terra/databases/checkv/checkv-db-v1.5.tar.gz"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  # fail hard
  set -euo pipefail

  # get version
  checkv -h | grep -Po "^CheckV [^:]+" | sed -e "s/CheckV //" | tee "VERSION"

  # extract CheckV DB
  tar -xzf ~{checkv_db}
  untarred_checkv_db=$(basename ~{checkv_db} .tar.gz)

  # run CheckV referencing the CheckV DB delineated by $CHECKVDB
  checkv end_to_end \
    ~{assembly} checkv_results/ \
    -d ${untarred_checkv_db} \
    -t ~{cpu} 

  # compile per-base results; NOTE: interpretation assumes a single viral genome
  python <<CODE
  import csv
  from collections import defaultdict

  # Read the quality summary file and extract the header
  data = defaultdict(list)
  with open('checkv_results/quality_summary.tsv', 'r') as infile:
    reader = csv.reader(infile, delimiter='\t')
    headers = next(reader)
    for row in reader:
      for i, value in enumerate(row):
        header = headers[i]
        if header in {"contig_length", "contamination", "completeness", "gene_count"}:
          # Convert to float
          value = float(value)
          data[header].append(value)
    
  # Summarize the data
  total_len = sum(data['contig_length'])
  total_genes = sum(data['gene_count'])
  # Contamination and completeness are reported as percents, so dividing by total reestablishes that
  per_base_contam = sum([v['contamination'] * v['contig_length'] for v in data.values()]) / total_len
  per_base_completeness = sum([v['completeness'] * v['contig_length'] for v in data.values()]) / total_len

  # Write the summary data for output exposure
  with open('TOTAL_GENES', 'w') as outfile:
    outfile.write(str(int(total_genes)))
  with open('PER_BASE_CONTAMINATION', 'w') as outfile:
    outfile.write(str(per_base_contam))
  with open('PER_BASE_COMPLETENESS', 'w') as outfile:
    outfile.write(str(per_base_completeness))

  CODE
  >>>
  output {
    String checkv_version = read_string("VERSION")
    Float per_base_contamination = read_float("PER_BASE_CONTAMINATION")
    Float weighted_contig_completeness = read_float("PER_BASE_COMPLETENESS")
    Int total_genes = read_int("TOTAL_GENES")
    File checkv_summary = "checkv_results/quality_summary.tsv"
    File checkv_contamination = "checkv_results/contamination.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}