version 1.0

task snpeff {
  input {
    String samplename
    String organism # organism name recorded for the custom SnpEff genome entry

    File reference_fasta # Reference full genome assembly FASTA
    File reference_gff # Reference GFF depicting annotated regions
    String? query_genes # comma-delimited list of strings
    File? bedfile
    File vcf

    String feature_qualifier = "product,locus_tag" # comma-delimited GFF feature qualifier(s) to use for comparison to query gene
    Boolean exact_match = false # use an exact match for qualifier mapping (always case-sensitive)
    Boolean ambiguous_contig = false # relate bedfile to GFF and FASTA ambiguous

    String docker = "us-docker.pkg.dev/general-theiagen/snpeff:5.4c"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  # variants are extracted only when a subset of genes/regions is requested; output
  # file suffixes reflect whether the extracted subset or the full VCF was annotated
  Boolean extract = defined(query_genes) || defined(bedfile)
  String vcf_scope = if extract then "extracted" else "full"
  command <<<
    # fail hard
    set -euo pipefail

    # obtain version
    snpeff -version 2>/dev/null | head -n 1 | cut -f 2 | tee VERSION

    # cap the JVM heap for every downstream snpeff call
    export JAVA_TOOL_OPTIONS="-Xmx~{memory}G"

    # SnpEff genome IDs become config keys, so restrict the samplename to safe characters
    genome_id=$(echo "~{samplename}" | sed 's/[^A-Za-z0-9._-]/_/g')
    data_dir="$(pwd)/snpeff_data"
    mkdir -p "${data_dir}/${genome_id}"

    # build a config from the container's default (retaining its codon tables and
    # settings), pointing data.dir at the local database directory and registering
    # the sample's genome under the organism name
    snpeff_dir=$(dirname "$(command -v snpeff)")
    sed "s|^data.dir *=.*|data.dir = ${data_dir}/|" "${snpeff_dir}/snpEff.config" > snpEff.config
    echo "${genome_id}.genome : ~{organism}" >> snpEff.config

    # strip any embedded ##FASTA section from the GFF; SnpEff reads the sequences
    # from the reference FASTA instead
    python3 <<CODE
    with open("~{reference_gff}") as fh, open("${data_dir}/${genome_id}/genes.gff", "w") as out:
        for line in fh:
            if line.strip().lower() in {'##fasta', '## fasta'}:
                break
            out.write(line)
    CODE

    # stage the reference FASTA where SnpEff expects it; zcat -f decompresses
    # gzipped input and passes plain text through unchanged
    zcat -f ~{reference_fasta} > "${data_dir}/${genome_id}/sequences.fa"

    # extract a sub-VCF only when a subset of genes/regions is requested, then
    # annotate that extracted subset instead of the full VCF
    if ~{extract}; then
      theiagene extract_variants \
        --vcf ~{vcf} \
        --reference_gff "${data_dir}/${genome_id}/genes.gff" \
        ~{'--query_genes "' + query_genes + '"'} \
        ~{"--bedfile " + bedfile} \
        ~{if exact_match then "--exact_match" else ""} \
        ~{if ambiguous_contig then "--ambiguous_contig" else ""} \
        --feature_qualifier "~{feature_qualifier}" \
        --output ~{samplename}_extracted.vcf
      snpeff_vcf="~{samplename}_extracted.vcf"
    else
      # by default, annotate the VCF passed directly to the task
      snpeff_vcf="~{vcf}"
    fi

    # count variant records (non-header, non-empty lines) in the VCF to annotate;
    # zcat -f decompresses bgzipped input and passes plain text through unchanged
    variant_count=$(zcat -f "${snpeff_vcf}" | grep -v "^#" | grep -c "[^[:space:]]" || true)

    # only build the database and annotate when the VCF holds at least one variant;
    # otherwise SnpEff has nothing to work on, so warn and skip
    if [ "${variant_count}" -gt 0 ]; then
      # build the SnpEff database from the GFF3 and reference *genome* FASTA
      snpeff build \
        -c snpEff.config \
        -gff3 \
        -noCheckCds \
        -noCheckProtein \
        -noLog \
        -v \
        "${genome_id}"

      # annotate the VCF; -ud 0 disables upstream/downstream annotations so only
      # variants overlapping features are annotated
      snpeff ann \
        -c snpEff.config \
        -nodownload \
        -noLog \
        -ud 0 \
        -stats ~{samplename}_~{vcf_scope}_snpeff_summary.html \
        "${genome_id}" \
        "${snpeff_vcf}" \
        > ~{samplename}_~{vcf_scope}_snpeff.vcf
    else
      echo "WARNING: no variants detected in ~{vcf_scope} VCF" >&2
    fi
  >>>
  output {
    String snpeff_version = read_string("VERSION")
    File? snpeff_annotated_vcf = "~{samplename}_~{vcf_scope}_snpeff.vcf"
    File? snpeff_summary_html = "~{samplename}_~{vcf_scope}_snpeff_summary.html"
    File? snpeff_genes_txt = "~{samplename}_~{vcf_scope}_snpeff_summary.genes.txt"
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
