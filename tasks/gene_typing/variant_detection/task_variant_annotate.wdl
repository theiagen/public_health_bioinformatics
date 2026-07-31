version 1.0

task variant_annotate {
  input {
    String samplename

    File reference_fasta 
    File reference_gff # Reference GFF depicting annotated regions 
    String? query_genes # comma-delimited list of strings
    File? bedfile
    File vcf
    
    String feature_qualifier = "product" # GFF feature qualifier to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)
    Boolean ambiguous_contig = false # relate bedfile to GFF and FASTA ambiguous

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ensembl-vep:116.0-dev"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    # fail hard
    set -euo pipefail

    # extract a sub-VCF only when a subset of genes/regions is requested, then
    # annotate that extracted subset instead of the full VCF
    if ~{if (defined(query_genes) || defined(bedfile)) then "true" else "false"}; then
      theiagene extract_variants \
        --vcf ~{vcf} \
        ~{if defined(query_genes) then "--query_genes ~{query_genes}" else ""} \
        ~{if defined(reference_gff) then "--reference_gff ~{reference_gff}" else ""} \
        ~{if defined(bedfile) then "--bedfile ~{bedfile}" else ""} \
        ~{if exact_match then "--exact_match" else ""} \
        ~{if ambiguous_contig then "--ambiguous_contig" else ""} \
        --feature_qualifier ~{feature_qualifier} \
        --output ~{samplename}.genes.vcf
      vep_vcf="~{samplename}.genes.vcf"
    else
      # by default, annotate the VCF passed directly to the task
      vep_vcf="~{vcf}"
    fi

# VEP requires the GFF sorted by contig then start coordinate; sort while
# keeping comment/header lines (leading '#') at the top of the file
python3 <<'PYEOF'
with open("~{reference_gff}") as fh:
    header = []
    records = []
    for line in fh:
        if line.startswith("#"):
            header.append(line)
        elif line.strip():
            records.append(line)
# sort by contig_id (col 1) then integer start coordinate (col 4)
records.sort(key=lambda l: (l.split("\t")[0], int(l.split("\t")[3])))
with open("reference_sorted.gff", "w") as out:
    out.writelines(header)
    out.writelines(records)
PYEOF

    # bgzip the sorted GFF (force overwrite if a stale .gz exists)
    bgzip -f reference_sorted.gff

    # index the bgzipped GFF with tabix
    tabix -p gff reference_sorted.gff.gz

    # bgzip the reference FASTA (to a local copy; inputs may be read-only) and
    # index with samtools faidx, which creates the .fai and .gzi indexes VEP needs
    # (FASTA files are indexed with faidx, not tabix)
    bgzip -c ~{reference_fasta} > ~{reference_fasta}.gz
    samtools faidx ~{reference_fasta}.gz

    vep \
      -i "${vep_vcf}" \
      --fasta ~{reference_fasta}.gz \
      --gff reference_sorted.gff.gz \
      --tab \
      --hgvs \
      --hgvsg \
      --hgvsp_use_prediction \
      --distance 0 \
      -o ~{samplename}_variant_annotations

    mv ~{samplename}_variant_annotations ~{samplename}_variant_annotations.tsv
    gunzip reference_sorted.gff.gz

    theiagene report_variants \
      --vcf "${vep_vcf}" \
      --vep_tsv ~{samplename}_variant_annotations.tsv \
      --reference_gff reference_sorted.gff \
      --suppress coding_sequence_variant,intergenic_variant \
      --feature_qualifier ~{feature_qualifier} \
      > VARIANT_REPORT
  >>>
  output {
    File variant_annotations_tsv = "~{samplename}_variant_annotations.tsv"
    File variant_annotation_warnings = "~{samplename}_variant_annotations_warnings.txt"
    File variant_annotations_html = "~{samplename}_variant_annotations_summary.html"
    File? variant_annotations_gene_vcf = "~{samplename}.genes.vcf"
    String variant_hits = read_string("VARIANT_REPORT")
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
