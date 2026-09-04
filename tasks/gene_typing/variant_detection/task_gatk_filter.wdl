version 1.0

task gatk_filter {
  input {
    String samplename
    File reference_genome
    File gvcf
    File gvcf_index

    # defaults informed by here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    Float min_variant_quality = 30
    Int min_depth = 10
    Float min_map_quality = 40
    Float min_quality_by_depth = 2.0
    Float max_fisher_strand_bias = 60
    Float max_strand_odds_ratio = 3
    String? filter_expression

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/gatk:4.6.2.0"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # fail hard
    set -euo pipefail

    # obtain version
    gatk --version | grep GATK | grep -Po "v[^ ]+$" | tee VERSION

    # create reference index and dictionary (have to symlink in case the source is not writeable)
    local_ref=$(basename ~{reference_genome})
    ln -s ~{reference_genome} ${local_ref}
    # index reference FASTA
    samtools faidx ${local_ref}
    gatk CreateSequenceDictionary -R ${local_ref}

    # assemble the VariantFiltration arguments, giving each filter its
    # own descriptive --filter-name
    python3 <<'CODE'
    # each entry is (filter_name, jexl_expression); VariantFiltration flags a
    # variant with filter_name when its expression evaluates to true. these
    # threshold expressions are written without spaces so bash word-splitting
    # keeps each token intact. the user-supplied filter_expression may contain
    # spaces, so it is handled separately as a quoted argument in the gatk call.
    # every threshold below has a task-level default, so all six always apply
    filters = [
        ('variant_quality_filter', 'QUAL<~{min_variant_quality}'),
        ('depth_filter', 'DP<~{min_depth}'),
        ('mapping_quality_filter', 'MQ<~{min_map_quality}'),
        ('quality_by_depth_filter', 'QD<~{min_quality_by_depth}'),
        ('fisher_strand_bias_filter', 'FS>~{max_fisher_strand_bias}'),
        ('strand_odds_ratio_filter', 'SOR>~{max_strand_odds_ratio}'),
    ]

    # build the flat argument list
    args = []
    for filter_name, filter_expr in filters:
        args += ["--filter-name", filter_name, "--filter-expression", filter_expr]

    with open("FILTER_EXPRESSION.txt", "w") as outfile:
        outfile.write(" ".join(args))

    print("Applying the following VariantFiltration filters:")
    for filter_name, filter_expr in filters:
        print('  {}: "{}"'.format(filter_name, filter_expr))

    # the user-supplied expression bypasses FILTER_EXPRESSION.txt, so report it here
    user_filter = '''~{default="" filter_expression}'''
    if user_filter:
        print('  user_filter: "{}"'.format(user_filter))
    CODE

    # call VariantFiltration with optional filter expression(s)
    mapfile -t FILTER_ARGS < FILTER_EXPRESSION.txt
    gatk --java-options "-Xmx~{memory}G" VariantFiltration \
      -R ${local_ref} \
      -V ~{gvcf} \
      -O ~{samplename}_filtered.g.vcf.gz \
      ~{'--filter-name "user_filter" --filter-expression "' + filter_expression + '"'} \
      ${FILTER_ARGS[@]}

    # call SelectVariants and drop those without PASS flags
    gatk --java-options "-Xmx~{memory}G" SelectVariants \
      -V ~{samplename}_filtered.g.vcf.gz \
      -O ~{samplename}_selected.g.vcf.gz \
      --exclude-filtered true

    # quantify the proportion of records that survived filtering; the FILTER
    # column (7) of every non-header record is either PASS or the name(s) of the
    # filter(s) it failed. records with no annotation (".") count toward the
    # denominator only. an empty VCF reports 0
    zcat ~{samplename}_filtered.g.vcf.gz \
      | awk -F'\t' '!/^#/ {total++; if ($7 == "PASS") pass++} END {printf "%.2f\n", (total > 0) ? (pass / total) * 100 : 0}' \
      | tee PERCENT_PASS.txt
  >>>
  output {
    String gatk_version = read_string("VERSION")
    File gatk_filtered_vcf = "~{samplename}_filtered.g.vcf.gz"
    File gatk_filtered_vcf_index = "~{samplename}_filtered.g.vcf.gz.tbi"
    File gatk_selected_vcf = "~{samplename}_selected.g.vcf.gz"
    File gatk_selected_vcf_index = "~{samplename}_selected.g.vcf.gz.tbi"
    Float gatk_percent_passing = read_float("PERCENT_PASS.txt")
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      preemptible: 0
      maxRetries: 3
  }
}
