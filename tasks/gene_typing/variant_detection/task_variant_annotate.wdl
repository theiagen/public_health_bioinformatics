version 1.0

task variant_annotate {
  input {
    String samplename

    File reference_gbff # Reference GBFF including annotated regions 
    String? query_genes # comma-delimited list of strings
    File vcf
    
    String feature_qualifier = "product" # GBFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/theiagene:1.0.0-dev"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # fail hard
    set -euo pipefail

    # annotate the protein-level consequences of variants overlapping the query genes;
    # this requires a reference GBFF (coding sequences), a VCF (variants) and query
    # genes, so it only runs when all three are provided. Kept non-fatal so an
    # annotation failure never discards the coverage outputs above.
    theiagene variant_annotation \
      --vcf ~{vcf} \
      --reference_gbff ~{reference_gbff} \
      ~{if defined(query_genes) then "--query_genes ~{query_genes}" else ""} \
      --feature_qualifier ~{feature_qualifier} \
      ~{if exact_match then "--exact_match" else ""} \
      --output ~{samplename}.variant_annotations.txt \
      || echo "WARNING: variant_annotation.py failed; continuing without a variant annotation report"

    if [ -f GENE_VARIANTS.vcf ]; then
      mv GENE_VARIANTS.vcf ~{samplename}.genes.vcf
    fi
  >>>
  output {
    File? gene_vcf = "~{samplename}.genes.vcf"
    String? variant_annotations = read_string("~{samplename}.variant_annotations.txt")
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
