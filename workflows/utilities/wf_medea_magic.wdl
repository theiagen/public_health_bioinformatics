version 1.0

import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as amr_search_task
import "../../tasks/alignment/task_bwa.wdl" as bwa_task
import "../../tasks/gene_typing/variant_detection/task_gatk_variants.wdl" as gatk_variants_task
import "../../tasks/gene_typing/variant_detection/task_gatk_filter.wdl" as gatk_filter_task
import "../../tasks/alignment/task_minimap2.wdl" as minimap2_task
import "../../tasks/gene_typing/variant_detection/task_clair3_variants.wdl" as clair3_task
import "../../tasks/utilities/data_handling/task_parse_mapping.wdl" as parse_mapping_task
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
    Float? cladetyper_max_distance
    # user-supplied reference fasta; when provided, overrides the hosted/organism reference
    File? reference_genome_fasta
    File? reference_gff
    # shared compute for the read_aligners (bwa for illumina, minimap2 for ont);
    String? read_aligner_docker
    Int? read_aligner_cpu
    Int? read_aligner_memory
    Int? read_aligner_disk_size
    # shared options for the variant callers (gatk variants/filter and clair3)
    String? gatk_docker
    Int? gatk_cpu
    Int? gatk_memory
    Int? gatk_disk_size
    # gatk-specific variant-calling options (illumina)
    Int? gatk_ploidy
    # gatk-specific filtering options (illumina)
    Float? gatk_filter_min_variant_quality
    Int? gatk_filter_min_depth
    Float? gatk_filter_min_map_quality
    Float? gatk_filter_min_quality_by_depth
    Float? gatk_filter_max_fisher_strand_bias
    Float? gatk_filter_max_strand_odds_ratio
    String? gatk_filter_expression
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
    File? query_genes_bed
    Boolean query_exact_match = false
    Int? min_gene_coverage_depth
    Int min_gene_coverage_map_quality = 40
  }
  # ORGANISM-SPECIFIC PARAMETER SETTING AND TASK CALLS
  if (medea_tag == "Candidozyma auris" || medea_tag == "Candida auris") {
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input:
        assembly_fasta = assembly,
        samplename = samplename,
        kmer_size = cladetyper_kmer_size,
        max_distance = cladetyper_max_distance,
        docker = cauris_cladetyper_docker_image
    }
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String cauris_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,B9J08_005340,B9J08_000401,B9J08_003102,B9J08_003737,B9J08_005343"
  }
  if (medea_tag == "Aspergillus fumigatus") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File afumigatus_variant_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fasta"
    # hosted GFF (user-supplied reference_gff takes precedence downstream) feeds gene coverage annotation
    File afumigatus_reference_gff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gff"
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String afumigatus_query_genes = "Cyp51A,HapE,AFUA_4G08340"
  }
  if (medea_tag == "Cryptococcus neoformans") {
    # hosted fasta (user-supplied fasta takes precedence downstream) feeds the variant-calling reference
    File cryptoneo_reference_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.fasta"
    # hosted GFF (user-supplied reference_gff takes precedence downstream) feeds gene coverage annotation
    File cryptoneo_reference_gff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gff"
    # organism-specific gene coverage targets used when query_genes is not user-supplied
    String cryptoneo_query_genes = "CNA00300"
  }

  # a user-supplied fasta takes precedence, otherwise the organism-specific reference is used
  # (cladetyper fasta for C. auris, hosted fasta for A. fumigatus and C. neoformans).
  # resolve the reference once; visible below (and in outputs) as File?
  String resolved_reference_fasta = select_first([reference_genome_fasta, cladetyper.assembly_reference, afumigatus_variant_fasta, cryptoneo_reference_fasta, ""])
  # Species-agnostic GFF resolution. A user-supplied reference_gff takes precedence,
  # otherwise the organism-specific reference is used (cladetyper outputs for
  # C. auris when a clade matches, hosted references for A. fumigatus and C. neoformans).
  # user cannot supply one reference file and not the other for GFF tasks to be functional
  if ((defined(reference_gff) && defined(reference_genome_fasta)) || (! defined(reference_gff) && ! defined(reference_genome_fasta))) {
    String resolved_reference_gff = select_first([reference_gff, cladetyper.annotated_reference_gff, afumigatus_reference_gff, cryptoneo_reference_gff, ""])
  }
  # The user-supplied query_genes takes priority; otherwise the
  # organism-specific default set (if any) is used. Inherently depends on variant calling
  String resolved_query_genes = select_first([query_genes, cauris_query_genes, afumigatus_query_genes, cryptoneo_query_genes, ""])


  # REFERENCE-BASED VARIANT CALLING
  # variant calling runs automatically whenever a reference fasta and read1 are available
  if (resolved_reference_fasta != "" && defined(read1)) {
    # Illumina short-read track: BWA alignment + GATK variant calling
    if (!ont_data) {
      call bwa_task.bwa as bwa_variant_calling {
        input:
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          reference_genome = resolved_reference_fasta,
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
          reference_genome = resolved_reference_fasta,
          ploidy = gatk_ploidy,
          docker = gatk_docker,
          cpu = gatk_cpu,
          memory = gatk_memory,
          disk_size = gatk_disk_size
      }
      call gatk_filter_task.gatk_filter as gatk_filter {
        input:
          samplename = samplename,
          reference_genome = resolved_reference_fasta,
          gvcf = gatk_variants.gatk_genotype_gvcf,
          gvcf_index = gatk_variants.gatk_genotype_gvcf_index,
          min_variant_quality = gatk_filter_min_variant_quality,
          min_depth = gatk_filter_min_depth,
          min_map_quality = gatk_filter_min_map_quality,
          min_quality_by_depth = gatk_filter_min_quality_by_depth,
          max_fisher_strand_bias = gatk_filter_max_fisher_strand_bias,
          max_strand_odds_ratio = gatk_filter_max_strand_odds_ratio,
          filter_expression = gatk_filter_expression,
          docker = gatk_docker,
          cpu = gatk_cpu,
          memory = gatk_memory,
          disk_size = gatk_disk_size
      }
    }
    # ONT long-read track: minimap2 alignment + Clair3 variant calling
    if (ont_data) {
      call minimap2_task.minimap2 as minimap2_variant_calling {
        input:
          query1 = select_first([read1]),
          reference = resolved_reference_fasta,
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
          reference_genome_file = resolved_reference_fasta,
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
    # Resolve newly-generated alignment files for downstream use
    File resolved_bam = select_first([bwa_variant_calling.sorted_bam, ont_bam_sorting.bam])
    File resolved_bai = select_first([bwa_variant_calling.sorted_bai, ont_bam_sorting.bai])
    File resolved_vcf = select_first([gatk_filter.gatk_selected_vcf, clair3_variant_calling.clair3_variants_vcf])
  }

  # GENE-CENTRIC COVERAGE CALCULATIONS AND VARIANT ANNOTATION
  # Terra can pass empty strings, so check for empty strings as well
  if ((select_first([resolved_reference_gff, ""]) != "" && resolved_query_genes != "") || defined(query_genes_bed)) {
    if (defined(resolved_bam)) {
      call gene_coverage_task.gene_coverage {
        input:
          bam = select_first([resolved_bam]),
          bai = select_first([resolved_bai]),
          samplename = samplename,
          reference_gff = resolved_reference_gff,
          bedfile = query_genes_bed,
          query_genes = resolved_query_genes,
          exact_match = query_exact_match,
          min_map_quality = min_gene_coverage_map_quality,
          min_depth = min_gene_coverage_depth
      }
      if (resolved_query_genes != "" && defined(resolved_vcf) && select_first([resolved_reference_gff, ""]) != "") {
        call variant_annotate_task.variant_annotate {
          input:
            samplename = samplename,
            reference_fasta = resolved_reference_fasta,
            reference_gff = select_first([resolved_reference_gff]),
            query_genes = resolved_query_genes,
            exact_match = query_exact_match,
            vcf = select_first([resolved_vcf]),
            bedfile = query_genes_bed
        }
      }
    }
  }

  # ANTIMICROBIAL RESISTANCE ANNOTATION
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
  # Output formatting
  if (resolved_reference_fasta != "") {
    File used_reference_fasta = resolved_reference_fasta
  }
  if (select_first([resolved_reference_gff, ""]) != "") {
    File used_reference_gff = select_first([resolved_reference_gff])
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
    # species-agnostic reference actually used (user-supplied or organism-derived)
    File? reference_fasta_used = used_reference_fasta
    File? reference_gff_used = used_reference_gff
    # variant calling - illumina (bwa alignment + gatk)
    String? bwa_version = bwa_variant_calling.bwa_version
    String? gatk_version = gatk_variants.gatk_version
    File? gatk_genotype_gvcf = gatk_variants.gatk_genotype_gvcf
    File? gatk_genotype_gvcf_index = gatk_variants.gatk_genotype_gvcf_index
    File? gatk_filtered_vcf = gatk_filter.gatk_filtered_vcf
    File? gatk_filtered_vcf_index = gatk_filter.gatk_filtered_vcf_index
    File? gatk_selected_vcf = gatk_filter.gatk_selected_vcf
    File? gatk_selected_vcf_index = gatk_filter.gatk_selected_vcf_index
    # variant calling - ont (minimap2 alignment + clair3)
    String? minimap2_version = minimap2_variant_calling.minimap2_version
    File? variant_calling_bam = resolved_bam
    File? variant_calling_bai = resolved_bai
    String? clair3_version = clair3_variant_calling.clair3_version
    File? clair3_variants_vcf = clair3_variant_calling.clair3_variants_vcf
    File? clair3_variants_gvcf = clair3_variant_calling.clair3_variants_gvcf
    String? clair3_variants_docker = clair3_variant_calling.clair3_variants_docker_image
    String? clair3_model_used = clair3_variant_calling.clair3_model_used
    # gene coverage
    File? gene_coverage_stats = gene_coverage.gene_coverage_stats
    String? gene_coverage_reads_mapped = gene_coverage.num_reads_aligned_to_genes
    String? gene_coverage_mean_breadth = gene_coverage.mean_breadth
    String? gene_coverage_mean_depth = gene_coverage.mean_depth
    Map[String, String]? gene_coverage_depth_by_gene = gene_coverage.depth_by_gene
    Map[String, String]? gene_coverage_breadth_by_gene = gene_coverage.breadth_by_gene
    Map[String, String]? gene_coverage_reads_by_gene = gene_coverage.reads_by_gene
    File? variant_annotation_warnings = variant_annotate.variant_annotation_warnings
    File? variant_annotation_summary = variant_annotate.variant_annotation_html
    File? variant_annotation_gene_vcf = variant_annotate.variant_annotation_gene_vcf
    File? variant_annotation_tsv = variant_annotate.variant_annotation_tsv
    String? variant_annotations = variant_annotate.variant_annotation
  }
}
