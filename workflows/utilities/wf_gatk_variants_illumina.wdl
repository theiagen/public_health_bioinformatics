version 1.0

import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/gene_typing/variant_detection/task_gatk_variants.wdl" as gatk_variants_task
import "../../tasks/gene_typing/variant_detection/task_gatk_filter.wdl" as gatk_filter_task

# GATK variant calling for Illumina paired-end data: BWA alignment, HaplotypeCaller
# genotyping, and hard filtering. Currently only used as a sub-workflow of medea_magic.
workflow gatk_variants_illumina {
  meta {
    description: "Align Illumina paired-end reads with BWA and call and filter variants with GATK"
  }
  input {
    File read1
    File read2 # paired-end reads are required
    File reference_genome_file
    String samplename
    # bwa runtime options
    String? bwa_docker
    Int? bwa_cpu
    Int? bwa_memory
    Int? bwa_disk_size
    # gatk variant-calling options
    Int? gatk_ploidy
    Boolean? gatk_decompress # output uncompressed genotype/filtered VCFs instead of bgzipped VCFs; the selected VCF is always uncompressed
    # gatk filtering options
    Float? gatk_filter_min_variant_quality
    Int? gatk_filter_min_depth
    Float? gatk_filter_min_map_quality
    Float? gatk_filter_min_quality_by_depth
    Float? gatk_filter_max_fisher_strand_bias
    Float? gatk_filter_max_strand_odds_ratio
    String? gatk_filter_expression
    # gatk runtime options, shared by the variant-calling and filtering tasks
    String? gatk_docker
    Int? gatk_cpu
    Int? gatk_memory
    Int? gatk_disk_size
  }
  call bwa_task.bwa {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      reference_genome = reference_genome_file,
      docker = bwa_docker,
      cpu = bwa_cpu,
      memory = bwa_memory,
      disk_size = bwa_disk_size
  }
  call gatk_variants_task.gatk_variants {
    input:
      samplename = samplename,
      bam = bwa.sorted_bam,
      bai = bwa.sorted_bai,
      reference_genome = reference_genome_file,
      ploidy = gatk_ploidy,
      decompress = gatk_decompress,
      docker = gatk_docker,
      cpu = gatk_cpu,
      memory = gatk_memory,
      disk_size = gatk_disk_size
  }
  call gatk_filter_task.gatk_filter {
    input:
      samplename = samplename,
      reference_genome = reference_genome_file,
      vcf = gatk_variants.gatk_genotype_vcf,
      vcf_index = gatk_variants.gatk_genotype_vcf_index,
      min_variant_quality = gatk_filter_min_variant_quality,
      min_depth = gatk_filter_min_depth,
      min_map_quality = gatk_filter_min_map_quality,
      min_quality_by_depth = gatk_filter_min_quality_by_depth,
      max_fisher_strand_bias = gatk_filter_max_fisher_strand_bias,
      max_strand_odds_ratio = gatk_filter_max_strand_odds_ratio,
      filter_expression = gatk_filter_expression,
      decompress = gatk_decompress,
      docker = gatk_docker,
      cpu = gatk_cpu,
      memory = gatk_memory,
      disk_size = gatk_disk_size
  }
  output {
    # read alignment
    String bwa_version = bwa.bwa_version
    File gatk_aligned_bam = bwa.sorted_bam
    File gatk_aligned_bai = bwa.sorted_bai
    # variant calling
    String gatk_version = gatk_variants.gatk_version
    File gatk_genotype_vcf = gatk_variants.gatk_genotype_vcf
    File gatk_genotype_vcf_index = gatk_variants.gatk_genotype_vcf_index
    # variant filtering
    File gatk_filtered_vcf = gatk_filter.gatk_filtered_vcf
    File gatk_filtered_vcf_index = gatk_filter.gatk_filtered_vcf_index
    File gatk_selected_vcf = gatk_filter.gatk_selected_vcf
    File gatk_selected_vcf_index = gatk_filter.gatk_selected_vcf_index
    Float gatk_percent_passing = gatk_filter.gatk_percent_passing
  }
}
