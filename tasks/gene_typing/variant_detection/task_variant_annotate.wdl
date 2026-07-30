version 1.0

task variant_annotate {
  input {
    String samplename

    File reference_fasta 
    File reference_gff # Reference GFF depicting annotated regions 
    String? query_genes # comma-delimited list of strings
    File vcf
    
    String feature_qualifier = "product" # GBFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ensembl-vep:116.0-dev"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # fail hard
    set -euo pipefail

    vep \
      -i ~{vcf} \
      --fasta ~{reference_fasta} \
      --gff ~{reference_gff} \
      --vcf \
      --hgvs \
      -o ~{samplename}_variant_annotations

    mv ~{samplename}_variant_annotations ~{samplename}_variant_annotations.tsv
  >>>
  output {
    File variant_annotations_tsv = "~{samplename}_variant_annotations.tsv"
    File variant_annotation_warnings = "~{samplename}_variant_annotations_warnings.txt"
    File variant_annotations_html = "~{samplename}_variant_annotations.html"
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
