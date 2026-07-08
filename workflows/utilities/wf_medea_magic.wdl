version 1.0

import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as amr_search_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/gene_typing/variant_detection/task_gatk_variants.wdl" as gatk_variants_task
import "../../tasks/gene_typing/variant_detection/task_gatk_filter.wdl" as gatk_filter_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task

workflow medea_magic {
  meta {
    description: "Workflow for fungal species typing and reference-based variant calling"
  }
  input {
    String samplename
    String medea_tag
    File assembly
    File? read1
    File? read2
    Boolean run_amr_search = false
    # variant calling logic
    Boolean run_variant_calling = true
    Boolean ont_data = false
    String? amr_search_docker_image
    String? cauris_cladetyper_docker_image
    # amr_search options
    Int? amr_search_cpu
    Int? amr_search_memory
    Int? amr_search_disk_size
    # cladetyper options - primarily files we host
    Int? cladetyper_kmer_size
    File? cladetyper_ref_clade1
    File? cladetyper_ref_clade1_annotated
    File? cladetyper_ref_clade2
    File? cladetyper_ref_clade2_annotated
    File? cladetyper_ref_clade3
    File? cladetyper_ref_clade3_annotated
    File? cladetyper_ref_clade4
    File? cladetyper_ref_clade4_annotated
    File? cladetyper_ref_clade5
    File? cladetyper_ref_clade5_annotated
    File? cladetyper_ref_clade6
    File? cladetyper_ref_clade6_annotated
    Float? cladetyper_max_distance
    # reference genomes - hosted fastas feed the variant-calling reference for each organism
    File afumigatus_reference_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fasta"
    File cryptoneo_reference_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.fasta"
    # user-supplied reference fasta; when provided, overrides the hosted/organism reference
    File? reference_genome_fasta
    # shared compute for the read_aligners (bwa for illumina, minimap2 for ont);
    # these tracks are mutually exclusive (ont_data), and each task still falls back
    # to its own docker default when read_aligner_docker is left unset
    String? read_aligner_docker
    Int? read_aligner_cpu
    Int? read_aligner_memory
    Int? read_aligner_disk_size
    # shared compute for the variant callers (gatk variants/filter and clair3);
    # these tracks are mutually exclusive (ont_data), and each task still falls back
    # to its own docker default when variant_caller_docker is left unset
    String? variant_caller_docker
    Int? variant_caller_cpu
    Int? variant_caller_memory
    Int? variant_caller_disk_size
    # gatk variant-calling options (illumina)
    Int? gatk_ploidy
    File? gatk_intervals_file
    # gatk filtering options (illumina)
    Int? gatk_filter_min_variant_quality
    Int? gatk_filter_min_depth
    Float? gatk_filter_min_map_quality
    Int? gatk_filter_min_quality_by_depth
    String? gatk_filter_expression
    # clair3 variant-calling options (ont)
    String? clair3_model
    Int? clair3_variant_quality
    Boolean? clair3_include_all_contigs
    Boolean? clair3_enable_haploid_precise
    Boolean? clair3_disable_phasing
    Boolean? clair3_enable_gvcf
    Boolean? clair3_enable_long_indel
  }
  if (medea_tag == "Candidozyma auris" || medea_tag == "Candida auris") {
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input:
        assembly_fasta = assembly,
        samplename = samplename,
        kmer_size = cladetyper_kmer_size,
        ref_clade1 = cladetyper_ref_clade1,
        ref_clade1_annotated = cladetyper_ref_clade1_annotated,
        ref_clade2 = cladetyper_ref_clade2,
        ref_clade2_annotated = cladetyper_ref_clade2_annotated,
        ref_clade3 = cladetyper_ref_clade3,
        ref_clade3_annotated = cladetyper_ref_clade3_annotated,
        ref_clade4 = cladetyper_ref_clade4,
        ref_clade4_annotated = cladetyper_ref_clade4_annotated,
        ref_clade5 = cladetyper_ref_clade5,
        ref_clade5_annotated = cladetyper_ref_clade5_annotated,
        ref_clade6 = cladetyper_ref_clade6,
        ref_clade6_annotated = cladetyper_ref_clade6_annotated,
        max_distance = cladetyper_max_distance,
        docker = cauris_cladetyper_docker_image
    }
    # only proceed if cladetyper retrieves an annotated_reference (e.g. non-functional for clade VI)
    if (cladetyper.annotated_reference != "None") {
      # cladetyper fasta output feeds the variant-calling alignment reference for C. auris
      File cauris_variant_fasta = cladetyper.assembly_reference
    }
  }
  if (medea_tag == "Aspergillus fumigatus") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File afumigatus_variant_fasta = afumigatus_reference_fasta
  }
  if (medea_tag == "Cryptococcus neoformans") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File cryptoneo_variant_fasta = cryptoneo_reference_fasta
  }
  # Reference-based variant calling. A single reference FASTA accounts for every organism:
  # a user-supplied fasta takes precedence, otherwise the organism-specific reference is used
  # (cladetyper fasta for C. auris, hosted fasta for A. fumigatus and C. neoformans).
  Array[File] variant_calling_reference_fastas = select_all([reference_genome_fasta, cauris_variant_fasta, afumigatus_variant_fasta, cryptoneo_variant_fasta])
  if (run_variant_calling && length(variant_calling_reference_fastas) > 0 && defined(read1)) {
    # Illumina short-read track: BWA alignment + GATK variant calling
    if (!ont_data) {
      call bwa_task.bwa as bwa_variant_calling {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          reference_genome = variant_calling_reference_fastas[0],
          cpu = read_aligner_cpu,
          memory = read_aligner_memory,
          disk_size = read_aligner_disk_size,
          docker = read_aligner_docker
      }
      call gatk_variants_task.gatk_variants as gatk_variants {
        input:
          samplename = samplename,
          bam = bwa_variant_calling.sorted_bam,
          bai = bwa_variant_calling.sorted_bai,
          reference_genome = variant_calling_reference_fastas[0],
          ploidy = gatk_ploidy,
          intervals_file = gatk_intervals_file,
          docker = variant_caller_docker,
          cpu = variant_caller_cpu,
          memory = variant_caller_memory,
          disk_size = variant_caller_disk_size
      }
      call gatk_filter_task.gatk_filter as gatk_filter {
        input:
          samplename = samplename,
          reference_genome = variant_calling_reference_fastas[0],
          gvcf = gatk_variants.gatk_genotype_gvcf,
          gvcf_index = gatk_variants.gatk_genotype_gvcf_index,
          min_variant_quality = gatk_filter_min_variant_quality,
          min_depth = gatk_filter_min_depth,
          min_map_quality = gatk_filter_min_map_quality,
          min_quality_by_depth = gatk_filter_min_quality_by_depth,
          filter_expression = gatk_filter_expression,
          docker = variant_caller_docker,
          cpu = variant_caller_cpu,
          memory = variant_caller_memory,
          disk_size = variant_caller_disk_size
      }
    }
    # ONT long-read track: minimap2 alignment + Clair3 variant calling
    if (ont_data) {
      call minimap2_task.minimap2 as minimap2_variant_calling {
        input:
          query1 = select_first([read1]),
          reference = variant_calling_reference_fastas[0],
          samplename = samplename,
          mode = "map-ont",
          output_sam = true,
          long_read_flags = true,
          docker = read_aligner_docker,
          cpu = read_aligner_cpu,
          memory = read_aligner_memory,
          disk_size = read_aligner_disk_size
      }
      # convert the minimap2 SAM to the sorted, indexed BAM Clair3 expects
      call parse_mapping_task.sam_to_sorted_bam as clair3_sorted_bam {
        input:
          sam = minimap2_variant_calling.minimap2_out,
          samplename = samplename
      }
      # index the reference FASTA; Clair3 requires the accompanying .fai
      call clair3_task.clair3_variants as clair3_variant_calling {
        input:
          alignment_bam_file = clair3_sorted_bam.bam,
          alignment_bam_file_index = clair3_sorted_bam.bai,
          reference_genome_file = variant_calling_reference_fastas[0],
          sequencing_platform = "ont",
          samplename = samplename,
          clair3_model = clair3_model,
          variant_quality = clair3_variant_quality,
          include_all_contigs = clair3_include_all_contigs,
          enable_haploid_precise = clair3_enable_haploid_precise,
          disable_phasing = clair3_disable_phasing,
          enable_gvcf = clair3_enable_gvcf,
          enable_long_indel = clair3_enable_long_indel,
          docker = variant_caller_docker,
          memory = variant_caller_memory,
          cpu = variant_caller_cpu,
          disk_size = variant_caller_disk_size
      }
    }
  }
  # Running AMR Search
  if (run_amr_search) {
    # Map containing the taxon tag reported by typing paired with it's taxon code for AMR search.
    Map[String, String] taxon_code = {
      "Candida auris" : "498019",
      "Candidozyma auris" : "498019"
    }
    # Checks for a match to the AMR_Search available taxon codes
    if (medea_tag == "Candida auris" || medea_tag == "Candidozyma auris") {
      call amr_search_task.amr_search {
        input:
          input_fasta = assembly,
          samplename = samplename,
          amr_search_database = taxon_code[medea_tag],
          cpu = amr_search_cpu,
          memory = amr_search_memory,
          disk_size = amr_search_disk_size,
          docker = amr_search_docker_image
      }
    }
  }
  output {
    # AMR_Search
    File? amr_search_results = amr_search.amr_search_json_output
    File? amr_results_csv = amr_search.amr_search_output_csv
    File? amr_results_pdf = amr_search.amr_search_output_pdf
    String? amr_search_all_resistances = amr_search.amr_search_all_resistances
    String? amr_search_associated_resistances = amr_search.amr_search_associated_resistances
    String? amr_search_docker = amr_search.amr_search_docker_image
    String? amr_search_version = amr_search.amr_search_version
    # c auris
    String? clade_type = cladetyper.gambit_cladetype
    String? cladetyper_version = cladetyper.gambit_version
    String? cladetyper_docker_image = cladetyper.gambit_cladetyper_docker_image
    String? cladetype_annotated_ref = cladetyper.annotated_reference
    # variant calling - illumina (bwa alignment + gatk)
    String? bwa_version = bwa_variant_calling.bwa_version
    File? variant_calling_bam = bwa_variant_calling.sorted_bam
    File? variant_calling_bai = bwa_variant_calling.sorted_bai
    String? gatk_version = gatk_variants.gatk_version
    File? gatk_genotype_gvcf = gatk_variants.gatk_genotype_gvcf
    File? gatk_genotype_gvcf_index = gatk_variants.gatk_genotype_gvcf_index
    File? gatk_filtered_vcf = gatk_filter.gatk_filtered_vcf
    File? gatk_selected_vcf = gatk_filter.gatk_selected_vcf
    # variant calling - ont (minimap2 alignment + clair3)
    String? minimap2_version = minimap2_variant_calling.minimap2_version
    File? ont_variant_calling_bam = clair3_sorted_bam.bam
    File? ont_variant_calling_bai = clair3_sorted_bam.bai
    String? clair3_version = clair3_variant_calling.clair3_version
    File? clair3_variants_vcf = clair3_variant_calling.clair3_variants_vcf
    File? clair3_variants_gvcf = clair3_variant_calling.clair3_variants_gvcf
    String? clair3_docker = clair3_variant_calling.clair3_variants_docker_image
    String? clair3_model_used = clair3_variant_calling.clair3_model_used
  }
}
