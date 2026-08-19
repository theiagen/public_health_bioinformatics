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
      --feature_type ~{feature_type} \
      --feature_qualifier ~{feature_qualifier} \
      --min_depth ~{min_depth} \
      --min_base_quality ~{min_base_quality} \
      --min_mapping_quality ~{min_map_quality} \
      ~{if defined(query_genes) then "--query_genes ~{query_genes}" else ""} \
      ~{if exact_match then "--exact_match" else ""} \
      ~{if defined(bedfile) then "--bedfile ~{bedfile}" else ""} \
      ~{if defined(reference_gff) then "--reference_gff ~{reference_gff}" else ""} \
      ~{if ambiguous_contig then "--ambiguous_contig" else ""}

    # rename files
    mv COVERAGE_STATS.tsv ~{samplename}.coverage_stats.tsv

    python3 <<CODE
    import json

    def load(path):
      with open(path, "r") as f:
        return {k.upper(): v for k, v in json.load(f).items()}

    def render_total(value):
      # every value is summed as a float, but a whole total -- a read count,
      # above all -- should read as 140 rather than 140.0
      if value != "" and float(value).is_integer():
        return str(int(value))
      return str(value)

    # bases quantified per gene -- the denominator its depth and breadth were
    # taken over, and the weight that turns a per-gene mean into a per-base one
    gene2lengths = load("LENGTHS_DICT.json")

    # a gene averaged on a "base" basis is weighted by its quantified length, so
    # the result is the mean over every base rather than the mean of per-gene
    # values; on a "query" basis every gene counts once no matter how long it is.
    # Depth and breadth are per-base quantities, so they average per base; reads
    # are counted per gene, so they average per query.
    for key, av_val in {"COVERAGE": "base", "DEPTH": "base", "READS": "query"}.items():
      data_dict = load(f"{key}_DICT.json")

      # theiagene reports a gene that resolved no coordinates as "NA": it was
      # never measured, which is not a measured zero, so it contributes to
      # neither the mean nor the total
      measured = {
        gene: float(value) for gene, value in data_dict.items()
        if value != "NA" and gene2lengths.get(gene, "NA") != "NA"
      }

      if av_val == "base":
        weights = {gene: float(gene2lengths[gene]) for gene in measured}
      elif av_val == "query":
        weights = {gene: 1.0 for gene in measured}
      else:
        raise ValueError(f"averaging basis for {key} must be 'base' or 'query', got {av_val}")

      # an empty weight total means nothing was measured; report that as blank
      # rather than as a number, matching how a gene with no coordinates reports
      denominator = sum(weights.values())
      if denominator:
        mean_data = sum(measured[g] * weights[g] for g in measured) / denominator
        total_data = sum(measured.values())
      else:
        mean_data = ""
        total_data = ""

      with open(f"MEAN_{key}", "w") as f:
        f.write(str(mean_data))
      with open(f"TOTAL_{key}", "w") as f:
        f.write(render_total(total_data))

      # deprecated outputs v4.2.0; typed Float, so an unmeasured S gene falls
      # back to 0.0 rather than emitting "NA"
      if "~{organism}".lower() == "sars-cov-2":
        sc2_s_gene_data = measured.get("S", 0.0)
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
    String gene_reads_aligned = read_string("TOTAL_READS")
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
