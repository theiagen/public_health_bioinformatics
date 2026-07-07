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

    # sequentially prepare variant filter expression
    FILTER_EXPRESSION="~{filter_expression}"
    ~{if defined(min_variant_quality) then "FILTER_EXPRESSION=$FILTER_EXPRESSION && QUAL < ~{min_variant_quality}" else ""}
    ~{if defined(min_depth) then "FILTER_EXPRESSION=$FILTER_EXPRESSION && DP < ~{min_depth}" else ""}
    ~{if defined(min_map_quality) then "FILTER_EXPRESSION=$FILTER_EXPRESSION && MQ < ~{min_map_quality}" else ""}
    ~{if defined(min_quality_by_depth) then "FILTER_EXPRESSION=$FILTER_EXPRESSION && QD < ~{min_quality_by_depth}" else ""}
    if [[ ! -z "$FILTER_EXPRESSION" ]]; then
      echo 'Filtering variants with the following expression: "${FILTER_EXPRESSION}"'
      FILTER_EXPRESSION="-filter ${FILTER_EXPRESSION}"
    fi

    # call VariantFiltration
    gatk --java-options "-Xmx~{memory}G" VariantFiltration \
      -R ~{reference_genome} \
      -V ~{samplename}_genotype.g.vcf.gz \
      -O ~{samplename}_filtered.g.vcf.gz \
      $FILTER_EXPRESSION

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
