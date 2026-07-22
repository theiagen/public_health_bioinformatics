version 1.0

import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as amr_search_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/gene_typing/variant_detection/task_freebayes.wdl" as freebayes_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
import "../../tasks/utilities/data_handling/task_fasta_utilities.wdl" as fasta_utilities_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/gene_typing/variant_detection/task_variant_annotate.wdl" as variant_annotate_task

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
    Boolean ont_data = false
    # amr_search options
    Int? amr_search_cpu
    Int? amr_search_memory
    Int? amr_search_disk_size
    String? amr_search_docker_image
    Boolean run_amr_search = false
    # cladetyper options - primarily files we host
    String? cauris_cladetyper_docker_image
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
    # user-supplied reference fasta; when provided, overrides the hosted/organism reference
    File? reference_genome_fasta
    File? reference_gbff
    # shared compute for the read_aligners (bwa for illumina, minimap2 for ont);
    String? read_aligner_docker
    Int? read_aligner_cpu
    Int? read_aligner_memory
    Int? read_aligner_disk_size
    # freebayes variant-calling & native input-filtering options (illumina)
    Int? freebayes_ploidy
    File? freebayes_targets_bed
    Int? freebayes_min_mapping_quality
    Int? freebayes_min_base_quality
    Float? freebayes_min_alternate_fraction
    Int? freebayes_min_alternate_count
    Int? freebayes_min_coverage
    Boolean? freebayes_enable_gvcf
    String? freebayes_docker
    Int? freebayes_cpu
    Int? freebayes_memory
    Int? freebayes_disk_size
    # clair3 variant-calling & filtering options (ont)
    String? clair3_model
    Int? clair3_variant_quality
    Boolean? clair3_include_all_contigs
    Boolean? clair3_enable_haploid_precise
    Boolean? clair3_disable_phasing
    Boolean? clair3_enable_gvcf
    Boolean? clair3_enable_long_indel
    String? clair3_docker
    Int? clair3_cpu
    Int? clair3_memory
    Int? clair3_disk_size
    # gene coverage options; user-supplied inputs take priority
    String? query_genes
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
    # cladetyper fasta output feeds the variant-calling alignment reference for C. auris
    File cauris_variant_fasta = cladetyper.assembly_reference
    # cladetyper GBFF output feeds gene coverage annotation, but only when a clade match is found
    if (cladetyper.annotated_reference != "None") {
      File cauris_variant_gbff = cladetyper.annotated_reference
    }
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String cauris_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,B9J08_005340,B9J08_000401,B9J08_003102,B9J08_003737,B9J08_005343"
  }
  if (medea_tag == "Aspergillus fumigatus") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File afumigatus_variant_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fasta"
    # hosted GBFF (user-supplied reference_gbff takes precedence downstream) feeds gene coverage annotation
    File afumigatus_variant_gbff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff"
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String afumigatus_query_genes = "Cyp51A,HapE,AFUA_4G08340"
  }
  if (medea_tag == "Cryptococcus neoformans") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File cryptoneo_variant_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.fasta"
    # hosted GBFF (user-supplied reference_gbff takes precedence downstream) feeds gene coverage annotation
    File cryptoneo_variant_gbff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gbff"
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String cryptoneo_query_genes = "CNA00300"
  }
  # Reference-based variant calling.
  # a user-supplied fasta takes precedence, otherwise the organism-specific reference is used
  # (cladetyper fasta for C. auris, hosted fasta for A. fumigatus and C. neoformans).
  Array[File] variant_calling_reference_fastas = select_all([reference_genome_fasta, cauris_variant_fasta, afumigatus_variant_fasta, cryptoneo_variant_fasta])
  # variant calling runs automatically whenever a reference fasta and read1 are available
  if (length(variant_calling_reference_fastas) > 0 && defined(read1)) {
    # Illumina short-read track: BWA alignment + freebayes variant calling
    File variant_calling_reference_fasta = variant_calling_reference_fastas[0]
    if (!ont_data) {
      call bwa_task.bwa as bwa_variant_calling {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          reference_genome = variant_calling_reference_fasta,
          cpu = read_aligner_cpu,
          memory = read_aligner_memory,
          disk_size = read_aligner_disk_size,
          docker = read_aligner_docker
      }
      # freebayes applies its native input filters (mapping/base quality, allele
      # support, coverage) during calling, so no separate post-call filter step is used
      call freebayes_task.freebayes as freebayes {
        input:
          samplename = samplename,
          bam = bwa_variant_calling.sorted_bam,
          bai = bwa_variant_calling.sorted_bai,
          reference_genome = variant_calling_reference_fasta,
          ploidy = freebayes_ploidy,
          targets_bed = freebayes_targets_bed,
          min_mapping_quality = freebayes_min_mapping_quality,
          min_base_quality = freebayes_min_base_quality,
          min_alternate_fraction = freebayes_min_alternate_fraction,
          min_alternate_count = freebayes_min_alternate_count,
          min_coverage = freebayes_min_coverage,
          output_gvcf = freebayes_enable_gvcf,
          docker = freebayes_docker,
          cpu = freebayes_cpu,
          memory = freebayes_memory,
          disk_size = freebayes_disk_size
      }
    }
    # ONT long-read track: minimap2 alignment + Clair3 variant calling
    if (ont_data) {
      call minimap2_task.minimap2 as minimap2_variant_calling {
        input:
          query1 = select_first([read1]),
          reference = variant_calling_reference_fasta,
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
      call parse_mapping_task.sam_to_sorted_bam as ont_bam_sorting {
        input:
          sam = minimap2_variant_calling.minimap2_out,
          samplename = samplename
      }
      # variant calling
      call clair3_task.clair3_variants as clair3_variant_calling {
        input:
          alignment_bam_file = ont_bam_sorting.bam,
          alignment_bam_file_index = ont_bam_sorting.bai,
          reference_genome_file = variant_calling_reference_fasta,
          sequencing_platform = "ont",
          samplename = samplename,
          clair3_model = clair3_model,
          variant_quality = clair3_variant_quality,
          include_all_contigs = clair3_include_all_contigs,
          enable_haploid_precise = clair3_enable_haploid_precise,
          disable_phasing = clair3_disable_phasing,
          enable_gvcf = clair3_enable_gvcf,
          enable_long_indel = clair3_enable_long_indel,
          docker = clair3_docker,
          memory = clair3_memory,
          cpu = clair3_cpu,
          disk_size = clair3_disk_size
      }
    }
  }
  # Gene coverage. The user-supplied query_genes takes priority; otherwise the
  # organism-specific default set (if any) is used. Inherently depends on variant calling
  Array[String] query_genes_options = select_all([query_genes, cauris_query_genes, afumigatus_query_genes, cryptoneo_query_genes])
  if (length(query_genes_options) > 0) {
    String resolved_query_genes = query_genes_options[0]
  }
  # a user-supplied reference_gbff takes precedence, otherwise the organism-specific GBFF
  # is used (cladetyper GBFF for C. auris when a clade matches, hosted GBFF for A. fumigatus and C. neoformans)
  Array[File] reference_gbff_options = select_all([reference_gbff, cauris_variant_gbff, afumigatus_variant_gbff, cryptoneo_variant_gbff])
  if (length(reference_gbff_options) > 0) {
    File resolved_reference_gbff = reference_gbff_options[0]
  }
  # tracks are mutually exclusive (ont_data), so each select_all yields at most one element
  Array[File] gene_coverage_bams = select_all([bwa_variant_calling.sorted_bam, ont_bam_sorting.bam])
  Array[File] gene_coverage_bais = select_all([bwa_variant_calling.sorted_bai, ont_bam_sorting.bai])
  Array[File] gene_coverage_vcfs = select_all([freebayes.freebayes_vcf, clair3_variant_calling.clair3_variants_vcf])
  if (length(gene_coverage_vcfs) > 0) {
    File gene_coverage_vcf = gene_coverage_vcfs[0]
  }
  if (length(gene_coverage_bams) > 0) {
    call gene_coverage_task.gene_coverage {
      input:
        bam = gene_coverage_bams[0],
        bai = gene_coverage_bais[0],
        samplename = samplename,
        reference_gbff = resolved_reference_gbff,
        query_genes = resolved_query_genes
    }
    if (defined(resolved_query_genes)) {
      call variant_annotate_task.variant_annotate {
        input:
          samplename = samplename,
          reference_gbff = select_first([resolved_reference_gbff]),
          query_genes = resolved_query_genes,
          vcf = select_first([gene_coverage_vcf])
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
    # variant calling - illumina (bwa alignment + freebayes)
    String? bwa_version = bwa_variant_calling.bwa_version
    File? variant_calling_bam = bwa_variant_calling.sorted_bam
    File? variant_calling_bai = bwa_variant_calling.sorted_bai
    String? freebayes_version = freebayes.freebayes_version
    String? freebayes_docker_image = freebayes.freebayes_docker_image
    File? freebayes_vcf = freebayes.freebayes_vcf
    File? freebayes_vcf_index = freebayes.freebayes_vcf_index
    File? freebayes_gvcf = freebayes.freebayes_gvcf
    File? freebayes_gvcf_index = freebayes.freebayes_gvcf_index
    # variant calling - ont (minimap2 alignment + clair3)
    String? minimap2_version = minimap2_variant_calling.minimap2_version
    File? ont_variant_calling_bam = ont_bam_sorting.bam
    File? ont_variant_calling_bai = ont_bam_sorting.bai
    String? clair3_version = clair3_variant_calling.clair3_version
    File? clair3_variants_vcf = clair3_variant_calling.clair3_variants_vcf
    File? clair3_variants_gvcf = clair3_variant_calling.clair3_variants_gvcf
    String? clair3_variants_docker = clair3_variant_calling.clair3_variants_docker_image
    String? clair3_model_used = clair3_variant_calling.clair3_model_used
    # gene coverage
    File? gene_coverage_stats = gene_coverage.gene_coverage_stats
    Map[String, Float]? gene_coverage_depth_by_gene = gene_coverage.depth_by_gene
    Map[String, Float]? gene_coverage_breadth_by_gene = gene_coverage.breadth_by_gene
    File? variant_annotation_gene_vcf = variant_annotate.gene_vcf
    String? variant_annotations = variant_annotate.variant_annotations
  }
}
