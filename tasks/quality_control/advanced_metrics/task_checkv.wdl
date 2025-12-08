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
  # fail hard, excluding pipefails
  set -eu

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
    -t ~{cpu} \
    || true

  if [ -e checkv_results/quality_summary.tsv ]; then
    echo "PASS" | tee CHECKV_STATUS
    echo "DEBUG: Extracting statistics"
    # col (2) is contig_length
    total_len=$(awk 'NR>1 {sum += $2} END {print sum}' checkv_results/quality_summary.tsv)
    # col (5) is gene_count
    awk 'NR>1 {sum += $5} END {print sum}' checkv_results/quality_summary.tsv | tee TOTAL_GENES

    # sum(col (2) contig_length * col (10) completeness) / total_len
    awk -v total_len="$total_len" 'NR>1 {sum += $2 * $10} END {result = sprintf("%.2f", sum / total_len); print result}' checkv_results/quality_summary.tsv | tee WEIGTHTED_COMPLETENESS

    # sum(col (2) contig_length * col (12) contamination) / total_len
    awk -v total_len="$total_len" 'NR>1 {sum += $2 * $12} END {result = sprintf("%.2f", sum / total_len); print result}' checkv_results/quality_summary.tsv | tee WEIGHTED_CONTAMINATION
  else
    echo "CheckV output files not detected"
    echo "FAIL" | tee CHECKV_STATUS
    touch WEIGHTED_CONTAMINATION
    touch WEIGHTED_COMPLETENESS
    touch TOTAL_GENES
  fi

  >>>
  output {
    String checkv_version = read_string("VERSION")
    String checkv_status = read_string("CHECKV_STATUS")
    String weighted_contamination = read_string("WEIGHTED_CONTAMINATION")
    String weighted_completeness = read_string("WEIGHTED_COMPLETENESS")
    String total_genes = read_string("TOTAL_GENES")
    File? checkv_summary = "checkv_results/quality_summary.tsv"
    File? checkv_contamination = "checkv_results/contamination.tsv"
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