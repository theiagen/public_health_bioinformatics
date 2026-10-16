version 1.0

task gene_coverage_theiagene {
  meta {
    description: "Gene depth and breadth of coverage with theiagene, behind the v4.3.0 gene_coverage interface"
  }
  input {
    File bam
    String samplename

    File? bai
    File snippy_outdir_tarball # Snippy writes reference/ref.gff from the GBFF it called variants against
    File? bedfile # BEDfile including region names and/or coordinates
    String? query_genes # comma-delimited list of strings

    String feature_type = "CDS" # GFF feature type to use for coordinate extraction
    String feature_qualifier = "product,locus_tag" # comma-delimited GFF attribute(s) to compare to query genes
    Boolean exact_match = false # use an exact match for qualifier mapping
    Boolean ambiguous_contig = false # apply coordinates from BED to first identified contig in BAM

    Int min_depth = 10 # minimum depth to count a base in breadth of coverage calculations
    Int min_quality = 0 # minimum base quality to count a base in breadth of coverage calculations

    # ensembl-vep:116.0 as published 2026-09-11 (theiagene 1.0.1, pysam 0.24.1); pinned by digest because
    # the 116.0 tag is rebuilt whenever theiagene is bumped
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ensembl-vep@sha256:d452689a44acf538f5471c300d16d5ad3fa67ddf8a5eb0d83191862bb19676cf"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    set -euo pipefail

    # Snippy's GFF shares the BAM's contig names (GBFF LOCUS), unlike NCBI GFFs (accession.version)
    tar -xzOf ~{snippy_outdir_tarball} --wildcards '*/reference/ref.gff' > reference.gff

    theiagene gene_coverage \
      --bam ~{bam} \
      --reference_gff reference.gff \
      --feature_type "~{feature_type}" \
      --feature_qualifier "~{feature_qualifier}" \
      --min_depth ~{min_depth} \
      --min_base_quality ~{min_quality} \
      ~{'--query_genes "' + query_genes + '"'} \
      ~{"--bedfile " + bedfile} \
      ~{if exact_match then "--exact_match" else ""} \
      ~{if ambiguous_contig then "--ambiguous_contig" else ""}

    # restore the v4.3.0 interface: resolved queries only in the maps (theiagene reports
    # unresolved queries as "NA"), and a three-column TSV reporting unresolved queries as 0
    python3 <<CODE
    import json

    for name in ("DEPTH_DICT.json", "COVERAGE_DICT.json"):
      with open(name, "r") as f:
        resolved = {k: float(v) for k, v in json.load(f).items() if k and v != "NA"}
      with open(name, "w") as f:
        # spoof Cromwell (Terra WDL) when nothing resolved, as v4.3.0 did
        json.dump(resolved or {"": 0}, f, indent=4)

    with open("COVERAGE_STATS.tsv", "r") as f:
      rows = [line.rstrip("\n").split("\t") for line in f if line.strip() and not line.startswith("#")]
    with open("~{samplename}.coverage_stats.tsv", "w") as out:
      out.write("#query\taverage_depth\tpercent_coverage")
      for query, depth, breadth, *_ in rows:
        out.write(f"\n{query}\t{0 if depth == 'NA' else depth}\t{0 if breadth == 'NA' else breadth}")
    CODE
  >>>
  output {
    File gene_coverage_stats = "~{samplename}.coverage_stats.tsv"
    Map[String, Float] depth_by_gene = read_json("DEPTH_DICT.json")
    Map[String, Float] breadth_by_gene = read_json("COVERAGE_DICT.json")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
