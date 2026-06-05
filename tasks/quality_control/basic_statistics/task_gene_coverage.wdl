version 1.0

task gene_coverage {
  input {
    File bam
    String samplename

    File? bai
    File? bedfile # BEDfile including region names and/or coordinates
    File? reference_gbff # GBFF including annotated regions 
    String? query_genes # comma-delimited list of strings
    
    String feature_type = "CDS" # GBFF feature type to use for coordinate extraction
    String feature_qualifier = "product" # GBFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)
    Boolean ambiguous_contig = false # apply coordinates from BED to first identified contig in BAM

    Int min_depth = 1 # minimum depth to count a base in breadth of coverage caclulations

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23.1-dev"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # fail hard
    set -euo pipefail

    # run calculations
    python3 /usr/bin/gene_coverage.py \
      --bam ~{bam} \
      --feature_type ~{feature_type} \
      --feature_qualifier ~{feature_qualifier} \
      --min_depth ~{min_depth} \
      ~{if defined(query_genes) then "--query_genes ~{query_genes}" else ""} \
      ~{if exact_match then "--exact_match" else ""} \
      ~{if defined(bedfile) then "--bedfile ~{bedfile}" else ""} \
      ~{if defined(reference_gbff) then "--reference_gbff ~{reference_gbff}" else ""} \
      ~{if ambiguous_contig then "--ambiguous_contig" else ""}

    mv COVERAGE_STATS.tsv ~{samplename}.coverage_stats.tsv
  >>>
  output {
    File gene_coverage_stats = "~{samplename}.coverage_stats.tsv"
    Map[String, Float] depth_by_gene = read_json("DEPTH_DICT.json")
    Map[String, Float] coverage_by_gene = read_json("COVERAGE_DICT.json")
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
