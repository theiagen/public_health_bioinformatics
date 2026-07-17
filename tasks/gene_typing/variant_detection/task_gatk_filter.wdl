version 1.0

task gatk_filter {
  input {
    String samplename
    File reference_genome
    File gvcf
    File gvcf_index

    Int? min_variant_quality
    Int? min_depth
    Float? min_map_quality
    Int? min_quality_by_depth
    String? filter_expression

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/gatk:4.6.2.0-dev"
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
    # spaces, so it is handled separately as a quoted argument in the gatk call
    filters = []

    # optional threshold-based filters, each with a name matching its function
    ~{if defined(min_variant_quality) then "filters.append(('variant_quality_filter', 'QUAL<~{min_variant_quality}'))" else ""}
    ~{if defined(min_depth) then "filters.append(('depth_filter', 'DP<~{min_depth}'))" else ""}
    ~{if defined(min_map_quality) then "filters.append(('mapping_quality_filter', 'MQ<~{min_map_quality}'))" else ""}
    ~{if defined(min_quality_by_depth) then "filters.append(('quality_by_depth_filter', 'QD<~{min_quality_by_depth}'))" else ""}

    # build the flat argument list; an empty list yields an empty file so no
    # filters are applied
    args = []
    for filter_name, filter_expr in filters:
        args += ["--filter-name", filter_name, "--filter-expression", filter_expr]

    with open("FILTER_EXPRESSION.txt", "w") as outfile:
        outfile.write(" ".join(args))

    if filters:
        print("Applying the following VariantFiltration filters:")
        for filter_name, filter_expr in filters:
            print('  {}: "{}"'.format(filter_name, filter_expr))
    else:
        print("No filter expression provided; running VariantFiltration without filters.")
    CODE

    # call VariantFiltration with optional filter expression(s)
    gatk --java-options "-Xmx~{memory}G" VariantFiltration \
      -R ${local_ref} \
      -V ~{gvcf} \
      -O ~{samplename}_filtered.g.vcf.gz \
      ~{if defined(filter_expression) then '--filter-name "user_filter" --filter-expression "~{filter_expression}"' else ""} \
      $(cat FILTER_EXPRESSION.txt)

    # call SelectVariants
    gatk --java-options "-Xmx~{memory}G" SelectVariants \
      -V ~{samplename}_filtered.g.vcf.gz \
      -O ~{samplename}_selected.g.vcf.gz
  >>>
  output {
    String gatk_version = read_string("VERSION")
    File gatk_filtered_vcf = "~{samplename}_filtered.g.vcf.gz"
    File gatk_selected_vcf = "~{samplename}_selected.g.vcf.gz"
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
