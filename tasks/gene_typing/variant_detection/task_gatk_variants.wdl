version 1.0

task gatk_variants {
  input {
    String samplename
    File bam
    File bai
    File reference_genome
    Int ploidy = 1 # integer indicating ploidy (N); default to haploid
    File? intervals_file

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

    # index reference FASTA
    samtools faidx ~{reference_genome}

    # create reference dictionary
    gatk CreateSequenceDictionary -R ~{reference_genome}

    # call HaplotypeCaller to generate an intermediate GVCF output depicting 
    # single-nucleotide polymorphisms (SNPs) and structural variants (SVs)
    # via local de-novo assembly of haplotypes in variant regions 
    gatk \
      --java-options "-Xms~{memory}G -Xmx~{memory}G" \
      HaplotypeCaller \
      -R ~{reference_genome} \
      -I ~{bam} \
      -O ~{samplename}_haplotypecall.g.vcf.gz \
      -ploidy ~{ploidy} \
      ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
      --tmp-dir . \
      -ERC GVCF

    # call GenotypeGVCFs
    # "This tool is able to handle any ploidy (or mix of ploidies) intelligently; 
    # there is no need to specify ploidy for non-diploid organisms."
    gatk --java-options "-Xmx~{memory}G" \
      GenotypeGVCFs \
      -R ~{reference_genome} \
      -V ~{samplename}_haplotypecall.g.vcf.gz \
      ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
      -O ~{samplename}_genotype.g.vcf.gz \

    # sequentially prepare variant filter expression and write it to FILTER_EXPRESSION.txt
    python3 <<'CODE'
    filter_conditions = []

    # start with the user-provided base filter expression, if any
    base_filter = "~{filter_expression}".strip()
    if base_filter:
        filter_conditions.append(base_filter)

    # append each optional threshold-based condition when its input is defined
    ~{if defined(min_variant_quality) then "filter_conditions.append('QUAL < ~{min_variant_quality}')" else ""}
    ~{if defined(min_depth) then "filter_conditions.append('DP < ~{min_depth}')" else ""}
    ~{if defined(min_map_quality) then "filter_conditions.append('MQ < ~{min_map_quality}')" else ""}
    ~{if defined(min_quality_by_depth) then "filter_conditions.append('QD < ~{min_quality_by_depth}')" else ""}

    # join all conditions with logical AND (empty string if none were provided)
    filter_expression = " && ".join(filter_conditions)

    with open("FILTER_EXPRESSION.txt", "w") as outfile:
        outfile.write(filter_expression)

    if filter_expression:
        print('Filtering variants with the following expression: "{}"'.format(filter_expression))
    else:
        print("No filter expression provided; running VariantFiltration without a filter.")
    CODE

    FILTER_EXPRESSION=$(cat FILTER_EXPRESSION.txt)

    # call VariantFiltration, applying the assembled filter only when one exists
    if [[ -n "$FILTER_EXPRESSION" ]]; then
      gatk --java-options "-Xmx~{memory}G" VariantFiltration \
        -R ~{reference_genome} \
        -V ~{samplename}_genotype.g.vcf.gz \
        -O ~{samplename}_filtered.g.vcf.gz \
        --filter-name "gatk_variant_filter" \
        -filter "$FILTER_EXPRESSION"
    else
      gatk --java-options "-Xmx~{memory}G" VariantFiltration \
        -R ~{reference_genome} \
        -V ~{samplename}_genotype.g.vcf.gz \
        -O ~{samplename}_filtered.g.vcf.gz
    fi

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
