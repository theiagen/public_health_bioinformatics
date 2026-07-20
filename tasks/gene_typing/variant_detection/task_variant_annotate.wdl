version 1.0

task variant_annotate {
  input {
    String samplename

    File? reference_gbff # Reference GBFF including annotated regions 
    File? reference_gff # Reference GFF including annotated regions
    File? reference_fa # Reference genome FASTA to be used in conjunction with GFF for gene coordinate extraction
    String query_genes # comma-delimited list of strings
    File vcf
    
    String feature_type = "CDS" # GBFF feature type to use for coordinate extraction
    String feature_qualifier = "product" # GBFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/pysam:1.23.2-dev"
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
    python3 /usr/bin/variant_annotation.py \
      --vcf ~{vcf} \
      ~{if defined(reference_gbff) then "--reference_gbff ~{reference_gbff}" else ""} \
      ~{if defined(reference_fa) then "--reference_fa ~{reference_fa}" else ""} \
      ~{if defined(reference_gff) then "--reference_gff ~{reference_gff}" else ""} \
      --query_genes ~{query_genes} \
      --feature_type ~{feature_type} \
      --feature_qualifier ~{feature_qualifier} \
      ~{if exact_match then "--exact_match" else ""} \
      --output ~{samplename}.variant_annotations.txt \
      || echo "WARNING: variant_annotation.py failed; continuing without a variant annotation report"
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
