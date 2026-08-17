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

    Int min_depth = 10 # minimum depth to count a base in breadth of coverage caclulations
    Int min_quality = 0 # minimum base quality to count a base in breadth of coverage caclulations

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
      --min_quality ~{min_quality} \
      ~{if defined(query_genes) then "--query_genes ~{query_genes}" else ""} \
      ~{if exact_match then "--exact_match" else ""} \
      ~{if defined(bedfile) then "--bedfile ~{bedfile}" else ""} \
      ~{if defined(reference_gff) then "--reference_gff ~{reference_gff}" else ""} \
      ~{if ambiguous_contig then "--ambiguous_contig" else ""} \

    # rename files
    mv COVERAGE_STATS.tsv ~{samplename}.coverage_stats.tsv

    python3 <<CODE
    import json
    for key in ["COVERAGE", "DEPTH", "READS"]:
      with open(f"{key}_DICT.json", "r") as f:
        data_dict = {k.upper(): v for k, v in json.load(f).items()}

      try:
        mean_data = sum(data_dict.values())/len(data_dict)
      except ZeroDivisionError:
        mean_data = ""

      with open(f"MEAN_{key}", "w") as f:
        f.write(mean_data)

      # deprecated outputs v4.2.0
      if "S" in data_dict and "~{organism}".lower() == "sars-cov-2":
        sc2_s_gene_data = data_dict["S"]
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
    String mean_reads = read_string("MEAN_READS")
    Map[String, Float] depth_by_gene = read_json("DEPTH_DICT.json")
    Map[String, Float] breadth_by_gene = read_json("COVERAGE_DICT.json")
    Map[String, Float] reads_by_gene = read_json("READS_DICT.json")
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
