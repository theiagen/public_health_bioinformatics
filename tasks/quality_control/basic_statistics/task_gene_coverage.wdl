version 1.0

task gene_coverage {
  input {
    File bam
    String samplename

    File? bai
    File? bedfile # BEDfile including region names and/or coordinates
    File? reference_gff # GFF including annotated regions
    String? query_genes # comma-delimited list of strings

    String feature_type = "CDS" # GFF feature type to use for coordinate extraction
    String feature_qualifier = "product" # GFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)
    Boolean ambiguous_contig = false # apply coordinates from BED to first identified contig in BAM

    Int min_depth = 10 # minimum depth to count a base
    Int min_map_quality = 40 # minimum mapping quality to count a base
    Int min_base_quality = 0 # minimum base quality to count a base

    String? organism # used to determine if S gene coverage should be reported for SARS-CoV-2

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ensembl-vep:116.0-dev"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # fail hard
    set -euo pipefail

    # run calculations
    theiagene gene_coverage \
      --bam ~{bam} \
      --feature_type "~{feature_type}" \
      --feature_qualifier "~{feature_qualifier}" \
      --min_depth ~{min_depth} \
      --min_base_quality ~{min_base_quality} \
      --min_mapping_quality ~{min_map_quality} \
      ~{'--query_genes "' + query_genes + '"'} \
      ~{if exact_match then "--exact_match" else ""} \
      ~{"--bedfile " + bedfile} \
      ~{"--reference_gff " + reference_gff} \
      ~{if ambiguous_contig then "--ambiguous_contig" else ""}

    # rename files
    mv COVERAGE_STATS.tsv ~{samplename}.coverage_stats.tsv

    # deprecated outputs v4.2.0; theiagene quantifies every gene, so the S gene
    # figures are read back out of its per-gene files rather than recomputed
    python3 <<CODE
    import json

    def load(path):
      with open(path, "r") as f:
        return {k.upper(): v for k, v in json.load(f).items()}

    # a gene theiagene reports as "NA" resolved to no coordinates and was never
    # measured; these outputs are typed Float, so an unmeasured or absent S gene
    # falls back to 0.0 rather than emitting "NA"
    for key in ("DEPTH", "COVERAGE"):
      s_gene_value = load(f"{key}_DICT.json").get("S", "NA")
      if "~{organism}".lower() == "sars-cov-2" and s_gene_value != "NA":
        sc2_s_gene_data = float(s_gene_value)
      else:
        sc2_s_gene_data = 0.0

      with open(f"SC2_S_GENE_{key}", "w") as f:
        f.write(str(sc2_s_gene_data))
    CODE
  >>>
  output {
    File gene_coverage_stats = "~{samplename}.coverage_stats.tsv"
    String mean_depth = read_string("MEAN_DEPTH")
    String mean_breadth = read_string("MEAN_COVERAGE")
    String num_reads_aligned_to_genes = read_string("TOTAL_READS")
    Map[String, String] depth_by_gene = read_json("DEPTH_DICT.json")
    Map[String, String] breadth_by_gene = read_json("COVERAGE_DICT.json")
    Map[String, String] reads_by_gene = read_json("READS_DICT.json")
    # deprecated v4.2.0
    Float sc2_s_gene_depth = read_string("SC2_S_GENE_DEPTH")
    Float sc2_s_gene_coverage = read_string("SC2_S_GENE_COVERAGE")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
