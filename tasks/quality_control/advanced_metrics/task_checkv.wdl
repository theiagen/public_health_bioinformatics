version 1.0

task checkv {
  meta {
    description: "Run CheckV on viral assemblies"
  }
  input {
    File assembly
    String samplename
    File checkv_db = "gs://theiagen-public-resources-rp/reference_data/databases/checkv/checkv-db-v1.5.tar.gz"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/checkv:1.0.3"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
  # fail hard
  set -euo pipefail

  # get version
  checkv -h | grep -Po "^CheckV [^:]+" | tee "VERSION"

  # extract CheckV DB
  echo "DEBUG: Extracting CheckV DB"
  tar -xzf ~{checkv_db}
  untarred_checkv_db=$(basename ~{checkv_db} .tar.gz)

  # run CheckV referencing the CheckV DB delineated by $CHECKVDB
  echo "DEBUG: Running CheckV end_to_end"
  checkv end_to_end \
    ~{assembly} checkv_results/ \
    -d ${untarred_checkv_db} \
    -t ~{cpu} 

  echo "DEBUG: Extracting statistics"
  # col (2) is contig_length
  total_len=$(awk 'NR>1 {sum += $2} END {print sum}' checkv_results/quality_summary.tsv)
  # col (5) is gene_count
  awk 'NR>1 {sum += $5} END {print sum}' checkv_results/quality_summary.tsv | tee TOTAL_GENES

  # sum(col (2) contig_length * col (10) completeness) / total_len
  awk -v total_len="$total_len" 'NR>1 {sum += $2 * $10} END {result = sprintf("%.2f", sum / total_len); print result}' checkv_results/quality_summary.tsv | tee WEIGTHTED_COMPLETENESS

  # sum(col (2) contig_length * col (12) contamination) / total_len
  awk -v total_len="$total_len" 'NR>1 {sum += $2 * $12} END {result = sprintf("%.2f", sum / total_len); print result}' checkv_results/quality_summary.tsv | tee WEIGHTED_CONTAMINATION

  >>>
  output {
    String checkv_version = read_string("VERSION")
    Float weighted_contamination = read_float("WEIGHTED_CONTAMINATION")
    Float weighted_completeness = read_float("WEIGTHTED_COMPLETENESS")
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