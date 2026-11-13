version 1.0

import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as amr_search_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task
import "../../tasks/quality_control/basic_statistics/task_mapping_stats.wdl" as mapping_stats_task
import "../../tasks/gene_typing/variant_detection/task_snpeff.wdl" as snpeff_task
import "../standalone_modules/wf_clair3_variants_ont.wdl" as clair3_variants_ont_workflow
import "wf_gatk_variants_illumina.wdl" as gatk_variants_illumina_workflow

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
    # for C. auris, supplying both reference_fasta and reference_gff consolidates the (by default bifurcated)
    # genome-wide and gene-level AMR variant calling paths onto the user-supplied reference
    File? reference_fasta
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
    Boolean? gatk_decompress # output uncompressed GATK genotype/filtered VCFs instead of bgzipped VCFs; the selected VCF is always uncompressed
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
    String? suppressed_variants # comma-delimited string of suppressed variants
    File? query_genes_bed
    Boolean query_exact_match = false
    Int min_gene_coverage_depth = 10 # mirror Snippy variants v4.6.0
    Int min_gene_coverage_map_quality = 60 # mirror Snippy variants v4.6.0
    Int min_gene_coverage_base_quality = 13 # mirror Snippy variants v4.6.0
  }
  # ORGANISM-SPECIFIC PARAMETER SETTING AND TASK CALLS
  if (medea_tag == "Candidozyma auris" || medea_tag == "Candida auris") {
    # C. auris reference selection is conducted via cladetyper
    call cauris_cladetyper.cauris_cladetyper as cladetyper {
      input:
        assembly_fasta = assembly,
        samplename = samplename,
        kmer_size = cladetyper_kmer_size,
        max_distance = cladetyper_max_distance,
        docker = cauris_cladetyper_docker_image
    }
    # B11205 (default C. auris AMR reference) defaults: product names, plus locus tags for genes annotated
    # only as "hypothetical protein" (MEC3, ERG3, PSK74852)
    String cauris_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase,FLO8,CJI82_01006,CJI82_00364,CJI82_03986"
  }
  if (medea_tag == "Aspergillus fumigatus") {
    File afumigatus_variant_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.fasta"
    File afumigatus_reference_gff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gff"
    String afumigatus_query_genes = "AFUA_4G06890,AFUA_6G05300,AFUA_4G08340"
  }
  if (medea_tag == "Cryptococcus neoformans") {
    File cryptoneo_reference_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.fasta"
    File cryptoneo_reference_gff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gff"
    String cryptoneo_query_genes = "sterol.14-demethylase"
  }

  # RESOLVE MULTI-SOURCE DOWNSTREAM INPUTS
  # a user-supplied fasta takes precedence, otherwise the organism-specific reference is used
  # (cladetyper fasta for C. auris, hosted fasta for A. fumigatus and C. neoformans).
  # resolve the reference once; visible below (and in outputs) as File?
  String resolved_reference_fasta = select_first([reference_fasta, cladetyper.assembly_reference, afumigatus_variant_fasta, cryptoneo_reference_fasta, ""])
  # Species-agnostic GFF resolution. A user-supplied reference_gff takes precedence,
  # otherwise the organism-specific reference is used (cladetyper outputs for
  # C. auris when a clade matches, hosted references for A. fumigatus and C. neoformans).
  # user cannot supply one reference file and not the other for GFF tasks to be functional
  if ((defined(reference_gff) && defined(reference_fasta)) || (! defined(reference_gff) && ! defined(reference_fasta))) {
    String resolved_reference_gff = select_first([reference_gff, cladetyper.annotated_reference_gff, afumigatus_reference_gff, cryptoneo_reference_gff, ""])
  }
  # The user-supplied query_genes takes priority; otherwise the organism-specific default set (if any) is used
  String resolved_query_genes = select_first([query_genes, cauris_query_genes, afumigatus_query_genes, cryptoneo_query_genes, ""])

  # REFERENCE-BASED VARIANT CALLING
  # variant calling runs automatically whenever a reference fasta and read1 are available;
  # Illumina data additionally requires read2, as only paired-end reads are supported
  if (resolved_reference_fasta != "" && defined(read1) && (ont_data || defined(read2))) {
    # Illumina paired-end track: BWA alignment + GATK variant calling (gatk_variants_illumina)
    if (!ont_data) {
      call gatk_variants_illumina_workflow.gatk_variants_illumina as gatk_genome_variants {
        input:
          read1 = select_first([read1]),
          read2 = select_first([read2]),
          reference_genome_file = resolved_reference_fasta,
          samplename = samplename,
          bwa_docker = read_aligner_docker,
          bwa_cpu = read_aligner_cpu,
          bwa_memory = read_aligner_memory,
          bwa_disk_size = read_aligner_disk_size,
          gatk_ploidy = gatk_ploidy,
          gatk_decompress = gatk_decompress,
          gatk_filter_min_variant_quality = gatk_filter_min_variant_quality,
          gatk_filter_min_depth = gatk_filter_min_depth,
          gatk_filter_min_map_quality = gatk_filter_min_map_quality,
          gatk_filter_min_quality_by_depth = gatk_filter_min_quality_by_depth,
          gatk_filter_max_fisher_strand_bias = gatk_filter_max_fisher_strand_bias,
          gatk_filter_max_strand_odds_ratio = gatk_filter_max_strand_odds_ratio,
          gatk_filter_expression = gatk_filter_expression,
          gatk_docker = gatk_docker,
          gatk_cpu = gatk_cpu,
          gatk_memory = gatk_memory,
          gatk_disk_size = gatk_disk_size
      }
    }
    # ONT long-read track: minimap2 alignment + Clair3 variant calling (Clair3_Variants_ONT)
    if (ont_data) {
      call clair3_variants_ont_workflow.clair3_variants_ont as clair3_genome_variants {
        input:
          read1 = select_first([read1]),
          reference_genome_file = resolved_reference_fasta,
          samplename = samplename,
          minimap2_docker = read_aligner_docker,
          minimap2_cpu = read_aligner_cpu,
          minimap2_memory = read_aligner_memory,
          minimap2_disk_size = read_aligner_disk_size,
          clair3_model = clair3_model,
          clair3_variant_quality = clair3_variant_quality,
          clair3_include_all_contigs = clair3_include_all_contigs,
          clair3_enable_haploid_precise = clair3_enable_haploid_precise,
          clair3_disable_phasing = clair3_disable_phasing,
          clair3_enable_gvcf = clair3_enable_gvcf,
          clair3_enable_long_indel = clair3_enable_long_indel,
          clair3_docker = clair3_docker,
          clair3_memory = clair3_memory,
          clair3_cpu = clair3_cpu,
          clair3_disk_size = clair3_disk_size
      }
    }
    # Resolve newly-generated alignment files for downstream use
    File resolved_bam = select_first([gatk_genome_variants.gatk_aligned_bam, clair3_genome_variants.clair3_aligned_bam])
    File resolved_bai = select_first([gatk_genome_variants.gatk_aligned_bai, clair3_genome_variants.clair3_aligned_bai])
    File resolved_vcf = select_first([gatk_genome_variants.gatk_selected_vcf, clair3_genome_variants.clair3_selected_vcf])
    # genome-wide read mapping statistics against the resolved reference
    call mapping_stats_task.mapping_stats {
      input:
        bamfile = resolved_bam,
        samplename = samplename,
        read1 = select_first([read1]),
        read2 = read2,
        vcf = resolved_vcf
    }
  }

  # C. AURIS AMR VARIANT CALLING
  # C. auris variant calling is bifurcated by default: genome-wide variants are called against the
  # cladetyper-selected (or user-supplied) reference above, while gene-level AMR characterization is standardized to B11205, so reads
  # are additionally called against it. These outputs are only used by the gene-centric modules; the AMR alignment
  # and its genome-wide variant files are not exposed.
  # Supplying both reference_fasta and reference_gff consolidates the two paths: the AMR variant calling is skipped
  # and the gene-centric modules use the genome-wide outputs and the user-supplied GFF. Supplying only one of the
  # two does not consolidate the paths, so gene-level AMR characterization remains standardized to B11205.
  if ((medea_tag == "Candidozyma auris" || medea_tag == "Candida auris") && !(defined(reference_fasta) && defined(reference_gff))) {
    if (defined(read1) && (ont_data || defined(read2))) {
      # bundle the AMR inputs so the gene-centric modules never mix references
      File cauris_amr_reference_fasta = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1-B11205_GCA_016772135.1_ASM1677213v1_genomic.fasta"
      File cauris_amr_reference_gff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1-B11205_GCA_016772135.1_ASM1677213v1_genomic.gff"
      if (!ont_data) {
        call gatk_variants_illumina_workflow.gatk_variants_illumina as cauris_amr_gatk_variants {
          input:
            read1 = select_first([read1]),
            read2 = select_first([read2]),
            reference_genome_file = cauris_amr_reference_fasta,
            samplename = samplename,
            bwa_docker = read_aligner_docker,
            bwa_cpu = read_aligner_cpu,
            bwa_memory = read_aligner_memory,
            bwa_disk_size = read_aligner_disk_size,
            gatk_ploidy = gatk_ploidy,
            gatk_decompress = gatk_decompress,
            gatk_filter_min_variant_quality = gatk_filter_min_variant_quality,
            gatk_filter_min_depth = gatk_filter_min_depth,
            gatk_filter_min_map_quality = gatk_filter_min_map_quality,
            gatk_filter_min_quality_by_depth = gatk_filter_min_quality_by_depth,
            gatk_filter_max_fisher_strand_bias = gatk_filter_max_fisher_strand_bias,
            gatk_filter_max_strand_odds_ratio = gatk_filter_max_strand_odds_ratio,
            gatk_filter_expression = gatk_filter_expression,
            gatk_docker = gatk_docker,
            gatk_cpu = gatk_cpu,
            gatk_memory = gatk_memory,
            gatk_disk_size = gatk_disk_size
        }
      }
      if (ont_data) {
        call clair3_variants_ont_workflow.clair3_variants_ont as cauris_amr_clair3_variants {
          input:
            read1 = select_first([read1]),
            reference_genome_file = cauris_amr_reference_fasta,
            samplename = samplename,
            minimap2_docker = read_aligner_docker,
            minimap2_cpu = read_aligner_cpu,
            minimap2_memory = read_aligner_memory,
            minimap2_disk_size = read_aligner_disk_size,
            clair3_model = clair3_model,
            clair3_variant_quality = clair3_variant_quality,
            clair3_include_all_contigs = clair3_include_all_contigs,
            clair3_enable_haploid_precise = clair3_enable_haploid_precise,
            clair3_disable_phasing = clair3_disable_phasing,
            clair3_enable_gvcf = clair3_enable_gvcf,
            clair3_enable_long_indel = clair3_enable_long_indel,
            clair3_docker = clair3_docker,
            clair3_memory = clair3_memory,
            clair3_cpu = clair3_cpu,
            clair3_disk_size = clair3_disk_size
        }
      }
      File cauris_amr_bam = select_first([cauris_amr_gatk_variants.gatk_aligned_bam, cauris_amr_clair3_variants.clair3_aligned_bam])
      File cauris_amr_bai = select_first([cauris_amr_gatk_variants.gatk_aligned_bai, cauris_amr_clair3_variants.clair3_aligned_bai])
      File cauris_amr_vcf = select_first([cauris_amr_gatk_variants.gatk_selected_vcf, cauris_amr_clair3_variants.clair3_selected_vcf])
    }
  }

  # GENE-CENTRIC INPUT RESOLUTION
  # the C. auris AMR variant calling outputs take priority when generated; otherwise (other organisms and
  # C. auris samples with a user-supplied reference_fasta and reference_gff) the genome-wide outputs are used
  String gene_reference_fasta = select_first([cauris_amr_reference_fasta, resolved_reference_fasta])
  String gene_reference_gff = select_first([cauris_amr_reference_gff, resolved_reference_gff, ""])
  if (gene_reference_gff != "") {
    File used_reference_gff = gene_reference_gff
  }

  # GENE-CENTRIC COVERAGE CALCULATIONS AND VARIANT ANNOTATION
  # Terra can pass empty strings, so check for empty strings as well
  if ((gene_reference_gff != "" && resolved_query_genes != "") || defined(query_genes_bed)) {
    if (defined(cauris_amr_bam) || defined(resolved_bam)) {
      call gene_coverage_task.gene_coverage {
        input:
          bam = select_first([cauris_amr_bam, resolved_bam]),
          bai = select_first([cauris_amr_bai, resolved_bai]),
          samplename = samplename,
          reference_gff = used_reference_gff,
          bedfile = query_genes_bed,
          query_genes = resolved_query_genes,
          exact_match = query_exact_match,
          min_map_quality = min_gene_coverage_map_quality,
          min_base_quality = min_gene_coverage_base_quality,
          min_depth = min_gene_coverage_depth,
          query_check_hard_fail = true
      }
      if (resolved_query_genes != "" && (defined(cauris_amr_vcf) || defined(resolved_vcf)) && gene_reference_gff != "") {
        call snpeff_task.snpeff as snpeff_gene_variants {
          input:
            samplename = samplename,
            organism = medea_tag,
            reference_fasta = gene_reference_fasta,
            reference_gff = gene_reference_gff,
            query_genes = resolved_query_genes,
            exact_match = query_exact_match,
            suppress = suppressed_variants,
            vcf = select_first([cauris_amr_vcf, resolved_vcf]),
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

  # OUTPUT FORMATTING
  if (resolved_reference_fasta != "") {
    File used_reference_fasta = resolved_reference_fasta
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
    # species-agnostic genome-wide reference actually used (user-supplied or organism-derived)
    File? reference_fasta_used = used_reference_fasta
    # GFF used by the gene-centric modules (the B11205 GFF for C. auris when its AMR variant calling runs)
    File? reference_gff_used = used_reference_gff
    # variant calling - illumina (bwa alignment + gatk)
    String? bwa_version = gatk_genome_variants.bwa_version
    String? gatk_genome_variants_version = gatk_genome_variants.gatk_version
    File? gatk_genome_variants_genotype_vcf = gatk_genome_variants.gatk_genotype_vcf
    File? gatk_genome_variants_genotype_vcf_index = gatk_genome_variants.gatk_genotype_vcf_index
    File? gatk_genome_variants_selected_vcf = gatk_genome_variants.gatk_selected_vcf
    File? gatk_genome_variants_selected_vcf_index = gatk_genome_variants.gatk_selected_vcf_index
    Float? gatk_genome_variants_percent_passing = gatk_genome_variants.gatk_percent_passing
    # variant calling - ont (minimap2 alignment + clair3)
    String? minimap2_version = clair3_genome_variants.minimap2_version
    File? variant_calling_bam = resolved_bam
    File? variant_calling_bai = resolved_bai
    String? clair3_genome_variants_version = clair3_genome_variants.clair3_version
    File? clair3_genome_variants_vcf = clair3_genome_variants.clair3_variants_vcf
    File? clair3_genome_variants_selected_vcf = clair3_genome_variants.clair3_selected_vcf
    File? clair3_genome_variants_gvcf = clair3_genome_variants.clair3_variants_gvcf
    String? clair3_genome_variants_docker = clair3_genome_variants.clair3_docker_image
    String? clair3_genome_variants_model_used = clair3_genome_variants.clair3_model_used
    # read mapping statistics
    String? mapping_stats_samtools_version = mapping_stats.samtools_version
    Float? mapping_stats_percentage_mapped_reads = mapping_stats.percentage_mapped_reads
    Float? mapping_stats_depth = mapping_stats.depth
    Float? mapping_stats_coverage = mapping_stats.coverage
    Int? mapping_stats_total_snps = mapping_stats.total_snps
    File? mapping_stats_stats = mapping_stats.stats
    File? mapping_stats_cov_hist = mapping_stats.cov_hist
    File? mapping_stats_cov_stats = mapping_stats.cov_stats
    File? mapping_stats_flagstat = mapping_stats.flagstat
    # gene query actually used (user-supplied or organism-derived)
    String gene_query = resolved_query_genes
    # gene coverage
    File? gene_coverage_stats = gene_coverage.gene_coverage_stats
    Map[String, String]? gene_coverage_depth_by_gene = gene_coverage.depth_by_gene
    Map[String, String]? gene_coverage_breadth_by_gene = gene_coverage.breadth_by_gene
    Map[String, String]? gene_coverage_reads_by_gene = gene_coverage.reads_by_gene
    String? gene_coverage_query_check = gene_coverage.gene_coverage_query_check
    String? gene_coverage_query_check_status = gene_coverage.gene_coverage_query_check_status
    # snpeff gene variant annotation
    String? snpeff_gene_variants_version = snpeff_gene_variants.snpeff_version
    File? snpeff_gene_variants_vcf = snpeff_gene_variants.snpeff_annotated_vcf
    File? snpeff_gene_variants_summary_html = snpeff_gene_variants.snpeff_summary_html
    File? snpeff_gene_variants_report_tsv = snpeff_gene_variants.snpeff_variant_report_tsv
    String? snpeff_gene_variants_hits = snpeff_gene_variants.snpeff_variant_hits
    String? snpeff_gene_variants_hits_abbreviated = snpeff_gene_variants.snpeff_variant_hits_abbreviated
  }
}
