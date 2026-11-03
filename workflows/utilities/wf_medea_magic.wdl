version 1.0

import "../../tasks/gene_typing/variant_detection/task_snippy_gene_query.wdl" as snippy_gene_query
import "../../tasks/gene_typing/variant_detection/task_snippy_variants.wdl" as snippy
import "../../tasks/species_typing/candidozyma/task_cauris_cladetyper.wdl" as cauris_cladetyper
import "../../tasks/gene_typing/drug_resistance/task_amr_search.wdl" as amr_search_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage_theiagene.wdl" as gene_coverage_theiagene_task

workflow medea_magic {
  meta {
    description: "Workflow for fungal species typing"
  }
  input {
    String samplename
    String medea_tag
    File assembly
    File? read1
    File? read2
    Boolean run_amr_search = false
    # subworkflow logic
    Boolean assembly_only = false
    String? amr_search_docker_image
    String? cauris_cladetyper_docker_image
    String? snippy_gene_query_docker_image
    String? snippy_variants_docker_image
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
    # snippy options - mostly files we host
    String? query_genes
    File afumigatus_reference_gbff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/aspergillus/Aspergillus_fumigatus_GCF_000002655.1_ASM265v1_genomic.gbff"
    File cryptoneo_reference_gbff = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/cryptococcus/Cryptococcus_neoformans_GCF_000091045.1_ASM9104v1_genomic.gbff"
    Int? snippy_map_qual
    Int? snippy_base_quality
    Int? snippy_min_coverage
    Float? snippy_min_frac
    Int? snippy_min_quality
    Int? snippy_maxsoft
    # gene coverage options
    String? gene_coverage_feature_type
    String? gene_coverage_feature_qualifier
    Int? gene_coverage_min_depth
    Int? gene_coverage_min_quality
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
    # default queries: product names where one gene resolves in every clade reference, otherwise
    # the matched clade reference's locus tags (FLO8, FLO8 N-terminus, MEC3, ERG3, PSK74852; from kzm-gatk-dev)
    String cauris_product_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase"
    if (cladetyper.gambit_cladetype == "Clade1") {
      String cauris_clade1_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase,B9J08_02359,B9J08_02360,B9J08_00960,B9J08_01595,B9J08_04226"
    }
    # Clade II writes products with "_" and annotates FLO8 as a single gene
    if (cladetyper.gambit_cladetype == "Clade2") {
      String cauris_clade2_query_genes = "FKS1,lanosterol_14-alpha_demethylase,uracil_phosphoribosyltransferase,sterol_24-C-methyltransferase,CJI96_0001121,CJI96_0001637,CJI96_0002270,CJI96_0004130"
    }
    if (cladetyper.gambit_cladetype == "Clade3") {
      String cauris_clade3_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase,CJI97_000401,CJI97_000402,CJI97_003175,CJI97_003811,CJI97_005426"
    }
    # Clade IV and V annotate FLO8 as a single gene
    if (cladetyper.gambit_cladetype == "Clade4") {
      String cauris_clade4_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase,CJJ07_003930,CJJ07_005307,CJJ07_003685,CJJ07_003405"
    }
    if (cladetyper.gambit_cladetype == "Clade5") {
      String cauris_clade5_query_genes = "FKS1,lanosterol.14-alpha.demethylase,uracil.phosphoribosyltransferase,sterol.24-C-methyltransferase,FDK38_002567,FDK38_000935,FDK38_001564,FDK38_004143"
    }
    # a user-supplied (e.g. Clade VI) annotation has unknown locus tags, so fall back to product names
    String cauris_query_genes = select_first([cauris_clade1_query_genes, cauris_clade2_query_genes, cauris_clade3_query_genes, cauris_clade4_query_genes, cauris_clade5_query_genes, cauris_product_query_genes])
    # only run snippy when cladetyper returns an annotation (none by default for Clade VI)
    if (defined(cladetyper.annotated_reference)) {
      if (!assembly_only) {
        call snippy.snippy_variants as snippy_cauris { # no ONT support right now
          input:
            reference_genome_file = select_first([cladetyper.annotated_reference]),
            read1 = select_first([read1]),
            read2 = read2,
            samplename = samplename,
            map_qual = snippy_map_qual,
            base_quality = snippy_base_quality,
            min_coverage = snippy_min_coverage,
            min_frac = snippy_min_frac,
            min_quality = snippy_min_quality,
            maxsoft = snippy_maxsoft,
            docker = snippy_variants_docker_image
        }
      }
      if (assembly_only) {
        call snippy.snippy_variants as snippy_cauris_assembly {
          input:
            reference_genome_file = select_first([cladetyper.annotated_reference]),
            assembly_fasta = assembly,
            samplename = samplename,
            map_qual = snippy_map_qual,
            base_quality = snippy_base_quality,
            min_coverage = snippy_min_coverage,
            min_frac = snippy_min_frac,
            min_quality = snippy_min_quality,
            maxsoft = snippy_maxsoft,
            docker = snippy_variants_docker_image
        }
      }
      call snippy_gene_query.snippy_gene_query as snippy_gene_query_cauris {
        input:
          samplename = samplename,
          snippy_variants_results = select_first([snippy_cauris.snippy_variants_results, snippy_cauris_assembly.snippy_variants_results]),
          reference = cladetyper.annotated_reference,
          query_gene = select_first([query_genes, cauris_query_genes]),
          docker = snippy_gene_query_docker_image
      }
      call gene_coverage_theiagene_task.gene_coverage_theiagene as gene_coverage_cauris {
        input:
          samplename = samplename,
          bam = select_first([snippy_cauris.snippy_variants_bam, snippy_cauris_assembly.snippy_variants_bam]),
          bai = select_first([snippy_cauris.snippy_variants_bai, snippy_cauris_assembly.snippy_variants_bai]),
          snippy_outdir_tarball = select_first([snippy_cauris.snippy_variants_outdir_tarball, snippy_cauris_assembly.snippy_variants_outdir_tarball]),
          query_genes = select_first([query_genes, cauris_query_genes]),
          feature_type = gene_coverage_feature_type,
          feature_qualifier = gene_coverage_feature_qualifier,
          min_depth = gene_coverage_min_depth,
          min_quality = gene_coverage_min_quality
      }
    }
  }
  if (medea_tag == "Aspergillus fumigatus") {
    # locus tags keep the search specific (from kzm-gatk-dev):
    # AFUA_4G06890 = cyp51A, AFUA_6G05300 = hapE, AFUA_4G08340 = cox10
    String afumigatus_query_genes = "AFUA_4G06890,AFUA_6G05300,AFUA_4G08340"
    if (!assembly_only) {
      call snippy.snippy_variants as snippy_afumigatus {
        input:
          reference_genome_file = afumigatus_reference_gbff,
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          map_qual = snippy_map_qual,
          base_quality = snippy_base_quality,
          min_coverage = snippy_min_coverage,
          min_frac = snippy_min_frac,
          min_quality = snippy_min_quality,
          maxsoft = snippy_maxsoft,
          docker = snippy_variants_docker_image
      }
    }
    if (assembly_only) {
      call snippy.snippy_variants as snippy_afumigatus_assembly {
        input:
          reference_genome_file = afumigatus_reference_gbff,
          assembly_fasta = assembly,
          samplename = samplename,
          map_qual = snippy_map_qual,
          base_quality = snippy_base_quality,
          min_coverage = snippy_min_coverage,
          min_frac = snippy_min_frac,
          min_quality = snippy_min_quality,
          maxsoft = snippy_maxsoft,
          docker = snippy_variants_docker_image
      }
    }
    call snippy_gene_query.snippy_gene_query as snippy_gene_query_afumigatus {
      input:
        samplename = samplename,
        snippy_variants_results = select_first([snippy_afumigatus.snippy_variants_results, snippy_afumigatus_assembly.snippy_variants_results]),
        reference = afumigatus_reference_gbff,
        query_gene = select_first([query_genes, afumigatus_query_genes]),
        docker = snippy_gene_query_docker_image
    }
    call gene_coverage_theiagene_task.gene_coverage_theiagene as gene_coverage_afumigatus {
      input:
        samplename = samplename,
        bam = select_first([snippy_afumigatus.snippy_variants_bam, snippy_afumigatus_assembly.snippy_variants_bam]),
        bai = select_first([snippy_afumigatus.snippy_variants_bai, snippy_afumigatus_assembly.snippy_variants_bai]),
        snippy_outdir_tarball = select_first([snippy_afumigatus.snippy_variants_outdir_tarball, snippy_afumigatus_assembly.snippy_variants_outdir_tarball]),
        query_genes = select_first([query_genes, afumigatus_query_genes]),
        feature_type = gene_coverage_feature_type,
        feature_qualifier = gene_coverage_feature_qualifier,
        min_depth = gene_coverage_min_depth,
        min_quality = gene_coverage_min_quality
    }
  }
  if (medea_tag == "Cryptococcus neoformans") {
    if (!assembly_only) {
      call snippy.snippy_variants as snippy_crypto {
        input:
          reference_genome_file = cryptoneo_reference_gbff,
          read1 = select_first([read1]),
          read2 = read2,
          samplename = samplename,
          map_qual = snippy_map_qual,
          base_quality = snippy_base_quality,
          min_coverage = snippy_min_coverage,
          min_frac = snippy_min_frac,
          min_quality = snippy_min_quality,
          maxsoft = snippy_maxsoft,
          docker = snippy_variants_docker_image
      }
    }
    if (assembly_only) {
      call snippy.snippy_variants as snippy_crypto_assembly {
        input:
          reference_genome_file = cryptoneo_reference_gbff,
          assembly_fasta = assembly,
          samplename = samplename,
          map_qual = snippy_map_qual,
          base_quality = snippy_base_quality,
          min_coverage = snippy_min_coverage,
          min_frac = snippy_min_frac,
          min_quality = snippy_min_quality,
          maxsoft = snippy_maxsoft,
          docker = snippy_variants_docker_image
      }
    }
    call snippy_gene_query.snippy_gene_query as snippy_gene_query_crypto {
      input:
        samplename = samplename,
        snippy_variants_results = select_first([snippy_crypto.snippy_variants_results, snippy_crypto_assembly.snippy_variants_results]),
        reference = cryptoneo_reference_gbff,
        query_gene = select_first([query_genes, "CNA00300"]), # CNA00300 is ERG11 for this reference genome
        docker = snippy_gene_query_docker_image
    }
    call gene_coverage_theiagene_task.gene_coverage_theiagene as gene_coverage_cryptoneo {
      input:
        samplename = samplename,
        bam = select_first([snippy_crypto.snippy_variants_bam, snippy_crypto_assembly.snippy_variants_bam]),
        bai = select_first([snippy_crypto.snippy_variants_bai, snippy_crypto_assembly.snippy_variants_bai]),
        snippy_outdir_tarball = select_first([snippy_crypto.snippy_variants_outdir_tarball, snippy_crypto_assembly.snippy_variants_outdir_tarball]),
        query_genes = select_first([query_genes, "CNA00300"]),
        feature_type = gene_coverage_feature_type,
        feature_qualifier = gene_coverage_feature_qualifier,
        min_depth = gene_coverage_min_depth,
        min_quality = gene_coverage_min_quality
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
    File? cladetype_annotated_ref = cladetyper.annotated_reference
    # snippy variants
    String snippy_variants_reference_genome = select_first([snippy_cauris.snippy_variants_reference_genome, snippy_cauris_assembly.snippy_variants_reference_genome, snippy_afumigatus.snippy_variants_reference_genome, snippy_afumigatus_assembly.snippy_variants_reference_genome, snippy_crypto.snippy_variants_reference_genome, snippy_crypto_assembly.snippy_variants_reference_genome, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_version = select_first([snippy_cauris.snippy_variants_version, snippy_cauris_assembly.snippy_variants_version,snippy_afumigatus.snippy_variants_version, snippy_afumigatus_assembly.snippy_variants_version, snippy_crypto.snippy_variants_version, snippy_crypto_assembly.snippy_variants_version, "No matching taxon detected"])
    String snippy_variants_query = select_first([snippy_gene_query_cauris.snippy_variants_query, snippy_gene_query_afumigatus.snippy_variants_query, snippy_gene_query_crypto.snippy_variants_query, "No matching taxon detected"])
    String snippy_variants_query_check = select_first([snippy_gene_query_cauris.snippy_variants_query_check, snippy_gene_query_afumigatus.snippy_variants_query_check, snippy_gene_query_crypto.snippy_variants_query_check, "No matching taxon detected"])
    String snippy_variants_hits = select_first([snippy_gene_query_cauris.snippy_variants_hits, snippy_gene_query_afumigatus.snippy_variants_hits, snippy_gene_query_crypto.snippy_variants_hits, "No matching taxon detected"])
    String snippy_variants_gene_query_results = select_first([snippy_gene_query_cauris.snippy_variants_gene_query_results, snippy_gene_query_afumigatus.snippy_variants_gene_query_results, snippy_gene_query_crypto.snippy_variants_gene_query_results, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_outdir_tarball = select_first([snippy_cauris.snippy_variants_outdir_tarball, snippy_cauris_assembly.snippy_variants_outdir_tarball, snippy_afumigatus.snippy_variants_outdir_tarball, snippy_afumigatus_assembly.snippy_variants_outdir_tarball, snippy_crypto.snippy_variants_outdir_tarball, snippy_crypto_assembly.snippy_variants_outdir_tarball, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_results = select_first([snippy_cauris.snippy_variants_results, snippy_cauris_assembly.snippy_variants_results,snippy_afumigatus.snippy_variants_results, snippy_afumigatus_assembly.snippy_variants_results, snippy_crypto.snippy_variants_results, snippy_crypto_assembly.snippy_variants_results, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_bam = select_first([snippy_cauris.snippy_variants_bam, snippy_cauris_assembly.snippy_variants_bam, snippy_afumigatus.snippy_variants_bam, snippy_afumigatus_assembly.snippy_variants_bam, snippy_crypto.snippy_variants_bam, snippy_crypto_assembly.snippy_variants_bam, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_bai = select_first([snippy_cauris.snippy_variants_bai, snippy_cauris_assembly.snippy_variants_bai, snippy_afumigatus.snippy_variants_bai, snippy_afumigatus_assembly.snippy_variants_bai, snippy_crypto.snippy_variants_bai, snippy_crypto_assembly.snippy_variants_bai, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_summary = select_first([snippy_cauris.snippy_variants_summary, snippy_cauris_assembly.snippy_variants_summary, snippy_afumigatus.snippy_variants_summary, snippy_afumigatus_assembly.snippy_variants_summary, snippy_crypto.snippy_variants_summary, snippy_crypto_assembly.snippy_variants_summary, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_num_reads_aligned = select_first([snippy_cauris.snippy_variants_num_reads_aligned, snippy_cauris_assembly.snippy_variants_num_reads_aligned, snippy_afumigatus.snippy_variants_num_reads_aligned, snippy_afumigatus_assembly.snippy_variants_num_reads_aligned, snippy_crypto.snippy_variants_num_reads_aligned, snippy_crypto_assembly.snippy_variants_num_reads_aligned, "No matching taxon detected"])
    String snippy_variants_coverage_tsv = select_first([snippy_cauris.snippy_variants_coverage_tsv, snippy_cauris_assembly.snippy_variants_coverage_tsv, snippy_afumigatus.snippy_variants_coverage_tsv, snippy_afumigatus_assembly.snippy_variants_coverage_tsv, snippy_crypto.snippy_variants_coverage_tsv, snippy_crypto_assembly.snippy_variants_coverage_tsv, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    String snippy_variants_num_variants = select_first([snippy_cauris.snippy_variants_num_variants, snippy_cauris_assembly.snippy_variants_num_variants, snippy_afumigatus.snippy_variants_num_variants, snippy_afumigatus_assembly.snippy_variants_num_variants, snippy_crypto.snippy_variants_num_variants, snippy_crypto_assembly.snippy_variants_num_variants, "No matching taxon detected"])
    String snippy_variants_percent_ref_coverage = select_first([snippy_cauris.snippy_variants_percent_ref_coverage, snippy_cauris_assembly.snippy_variants_percent_ref_coverage, snippy_afumigatus.snippy_variants_percent_ref_coverage, snippy_afumigatus_assembly.snippy_variants_percent_ref_coverage, snippy_crypto.snippy_variants_percent_ref_coverage, snippy_crypto_assembly.snippy_variants_percent_ref_coverage, "No matching taxon detected"])
    # gene coverage
    File gene_coverage_stats = select_first([gene_coverage_cauris.gene_coverage_stats, gene_coverage_afumigatus.gene_coverage_stats, gene_coverage_cryptoneo.gene_coverage_stats, "gs://theiagen-public-resources-rp/empty_files/no_match_detected.txt"])
    Map[String, Float] depth_by_gene = select_first([gene_coverage_cauris.depth_by_gene, gene_coverage_afumigatus.depth_by_gene, gene_coverage_cryptoneo.depth_by_gene, {"": 0}])
    Map[String, Float] percent_coverage_by_gene = select_first([gene_coverage_cauris.breadth_by_gene, gene_coverage_afumigatus.breadth_by_gene, gene_coverage_cryptoneo.breadth_by_gene, {"": 0}])
  }
}