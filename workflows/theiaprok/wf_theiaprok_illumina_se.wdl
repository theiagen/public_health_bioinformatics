version 1.0

import "../../workflows/utilities/wf_digger_denovo.wdl" as digger_denovo
import "../../tasks/gene_typing/annotation/task_bakta.wdl" as bakta_task
import "../../tasks/gene_typing/annotation/task_prokka.wdl" as prokka_task
import "../../tasks/gene_typing/drug_resistance/task_amrfinderplus.wdl" as amrfinderplus
import "../../tasks/gene_typing/drug_resistance/task_resfinder.wdl" as resfinder
import "../../tasks/gene_typing/plasmid_detection/task_plasmidfinder.wdl" as plasmidfinder_task
import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate_task
import "../../tasks/quality_control/advanced_metrics/task_busco.wdl" as busco_task
import "../../tasks/quality_control/advanced_metrics/task_mummer_ani.wdl" as ani_task
import "../../tasks/quality_control/basic_statistics/task_cg_pipeline.wdl" as cg_pipeline
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/species_typing/multi/task_ts_mlst.wdl" as ts_mlst_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/contamination/task_kmerfinder.wdl" as kmerfinder_task
import "../../tasks/taxon_id/task_gambit.wdl" as gambit_task
import "../../tasks/gene_typing/drug_resistance/task_gamma.wdl" as gamma_task
import "../../tasks/utilities/data_export/task_export_taxon_table.wdl" as export_taxon_table_task
import "../../tasks/utilities/data_handling/task_arln_stats.wdl" as arln_stats
import "../utilities/file_handling/wf_concatenate_illumina_lanes.wdl" as concatenate_lanes_workflow
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../utilities/wf_read_QC_trim_se.wdl" as read_qc

workflow theiaprok_illumina_se {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of paired-end bacterial NGS data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1

    # optional additional lanes
    File? read1_lane2
    File? read1_lane3
    File? read1_lane4

    Int? genome_length
    # export taxon table parameters
    String? run_id
    String? collection_date
    String? originating_lab
    String? city
    String? county
    String? zip
    File? taxon_tables
    String terra_project="NA"
    String terra_workspace="NA"
    # read screen paramaters
    Boolean skip_screen = false 
    Boolean skip_mash = false
    Int min_reads = 7472
    Int min_basepairs = 2241820
    Int min_genome_length = 100000
    Int max_genome_length = 18040666
    Int min_coverage = 10
    # trimming parameters
    Int trim_min_length = 25
    Int trim_quality_min_score = 30
    Int trim_window_size = 4
    # module options
    Boolean perform_characterization = true # by default run all characterization steps
    Boolean amrfinder_use_gff = false # by default use nucleotide fasta for amrfinderplus, but user can set this to true if they want to use a gff and protein fasta file
    Boolean call_ani = false # by default do not call ANI task, but user has ability to enable this task if working with enteric pathogens or supply their own high-quality reference genome
    Boolean call_kmerfinder = false 
    Boolean call_resfinder = false
    Boolean call_plasmidfinder = true
    Boolean call_abricate = false
    Boolean call_gamma = false
    Boolean call_arln_stats = false
    Boolean mlst_scheme_override = false # If true, will force E. coli scheme to be used when Gambit predicts Escherichia coli, otherwise will return scheme MLST predicts.
    Boolean mlst_run_secondary_scheme = false # If true, will run secondary scheme if primary scheme is of ecoli or abaumannii, these two have multiple schemes that are relevant.
    String abricate_db = "vfdb"
    String genome_annotation = "prokka" # options: "prokka" or "bakta"
    String bakta_db = "full" # Default: "light" or "full"
    String? expected_taxon # allow user to provide organism (e.g. "Clostridioides_difficile") string to amrfinder. Useful when gambit does not predict the correct species
    # qc check parameters
    File? qc_check_table
  }
  call versioning.version_capture {
    input:
  }
  if (defined(read1_lane2)) {
    call concatenate_lanes_workflow.concatenate_illumina_lanes {
      input:
        samplename = samplename,
        read1_lane1 = read1,
        read1_lane2 = select_first([read1_lane2]),
        read1_lane3 = read1_lane3,
        read1_lane4 = read1_lane4
    }
  }
  if (! skip_screen) {
    call screen.check_reads_se as raw_check_reads { 
      input:
        read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        skip_mash = skip_mash,
        expected_genome_length = genome_length,
        workflow_series = "theiaprok"
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    call read_qc.read_QC_trim_se as read_QC_trim {
      input:
        samplename = samplename,
        read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        workflow_series = "theiaprok"
    }
    if (! skip_screen) {
      call screen.check_reads_se as clean_check_reads {
        input:
          read1 = read_QC_trim.read1_clean,
          min_reads = min_reads,
          min_basepairs = min_basepairs,
          min_genome_length = min_genome_length,
          max_genome_length = max_genome_length,
          min_coverage = min_coverage,
          skip_mash = skip_mash,
          expected_genome_length = genome_length,
          workflow_series = "theiaprok"
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      call digger_denovo.digger_denovo {
        input:
          samplename = samplename,
          read1 = read_QC_trim.read1_clean
      }
      call quast_task.quast {
        input:
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename
      }
      call cg_pipeline.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
          samplename = samplename,
          genome_length = select_first([genome_length, quast.genome_length])
      }
      call cg_pipeline.cg_pipeline as cg_pipeline_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename,
          genome_length = select_first([genome_length, quast.genome_length])
      }
      call busco_task.busco {
          input:
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename
      }
      if (perform_characterization) {
        call gambit_task.gambit {
          input:
            assembly = digger_denovo.assembly_fasta,
            samplename = samplename
        }      
        if (call_ani) {
          call ani_task.animummer as ani {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename
          }
        }
        if (call_kmerfinder) {
          call kmerfinder_task.kmerfinder_bacteria as kmerfinder {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename
          }
        }
        call amrfinderplus.amrfinderplus_nuc as amrfinderplus_task {
          input:
            assembly = digger_denovo.assembly_fasta,
            annotation_assembly = select_first([prokka.prokka_fna,bakta.bakta_fna]),
            samplename = samplename,
            protein_fasta = select_first([prokka.prokka_faa,bakta.bakta_faa]),
            gff = select_first([prokka.prokka_gff,bakta.bakta_gff3]),
            organism = select_first([expected_taxon, gambit.gambit_predicted_taxon]),
            annotation_format = genome_annotation,
            use_gff = amrfinder_use_gff
        }
        if (call_gamma){
          call gamma_task.gamma{
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename
          }
        }
        if (call_resfinder) {
          call resfinder.resfinder as resfinder_task {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename,
              organism = select_first([expected_taxon, gambit.gambit_predicted_taxon])
          }
        }      
        call ts_mlst_task.ts_mlst {
          input: 
            assembly = digger_denovo.assembly_fasta,
            samplename = samplename,
            taxonomy = select_first([expected_taxon, gambit.gambit_predicted_taxon]),
            run_secondary_scheme = mlst_run_secondary_scheme,
            scheme_override = mlst_scheme_override
        }
        if (genome_annotation == "prokka") {
          call prokka_task.prokka {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename
          }
        }
        if (genome_annotation == "bakta") {  
          if (bakta_db == "light") {  
            File bakta_db_light = "gs://theiagen-public-resources-rp/reference_data/databases/bakta/bakta_db_light_2025-01-23.tar.gz"  
          }  
          if (bakta_db == "full") {  
            File bakta_db_full = "gs://theiagen-public-resources-rp/reference_data/databases/bakta/bakta_db_full_2024-01-23.tar.gz"            
          }  
          if (!(bakta_db == "light" || bakta_db == "full")) {  
              File bakta_custom_db = bakta_db  
          } 
          call bakta_task.bakta {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename,
              bakta_db_selected = select_first([bakta_custom_db, bakta_db_light, bakta_db_full])
          }
        }
        if (call_plasmidfinder) {
          call plasmidfinder_task.plasmidfinder {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename
          }
        }
        if (call_abricate) {
          call abricate_task.abricate {
            input:
              assembly = digger_denovo.assembly_fasta,
              samplename = samplename,
              database = abricate_db
          }
        }
        if (defined(qc_check_table)) {
          call qc_check.qc_check_phb as qc_check_task {
            input:
              qc_check_table = qc_check_table,
              expected_taxon = expected_taxon,
              gambit_predicted_taxon = gambit.gambit_predicted_taxon,
              qc_check_inputs = {
                "num_reads_raw1": select_first([read_QC_trim.fastq_scan_raw1, read_QC_trim.fastqc_raw1]),
                "num_reads_clean1": select_first([read_QC_trim.fastq_scan_clean1, read_QC_trim.fastqc_clean1]),
                "r1_mean_q_raw": cg_pipeline_raw.r1_mean_q,
                "r1_mean_readlength_raw": cg_pipeline_raw.r1_mean_readlength,
                "r1_mean_q_clean": cg_pipeline_clean.r1_mean_q,
                "r1_mean_readlength_clean": cg_pipeline_clean.r1_mean_readlength,   
                "est_coverage_raw": cg_pipeline_raw.est_coverage,
                "est_coverage_clean": cg_pipeline_clean.est_coverage,
                "midas_secondary_genus_abundance": read_QC_trim.midas_secondary_genus_abundance,
                "midas_secondary_genus_coverage": read_QC_trim.midas_secondary_genus_coverage,
                "assembly_length": quast.genome_length,
                "number_contigs": quast.number_contigs,
                "n50_value": quast.n50_value,
                "quast_gc_percent": quast.gc_percent,
                "busco_completeness": busco.busco_results,
                "ani_highest_percent": ani.ani_highest_percent,
                "ani_highest_percent_bases_aligned": ani.ani_highest_percent_bases_aligned
              }
          }
        }
        call merlin_magic_workflow.merlin_magic {
          input:
            merlin_tag = select_first([expected_taxon, gambit.merlin_tag]),
            assembly = digger_denovo.assembly_fasta,
            samplename = samplename,
            read1 = read_QC_trim.read1_clean,
            paired_end = false
        }
        if (defined(taxon_tables)) {
          call export_taxon_table_task.export_taxon_table {
            input:
              terra_project = terra_project,
              terra_workspace = terra_workspace,
              gambit_predicted_taxon = gambit.gambit_predicted_taxon,
              taxon_table = taxon_tables,
              samplename = samplename,
              columns_to_export = {
                "abricate_abaum_database": merlin_magic.abricate_abaum_database,
                "abricate_abaum_docker": merlin_magic.abricate_abaum_docker,
                "abricate_abaum_plasmid_tsv": merlin_magic.abricate_abaum_results,
                "abricate_abaum_plasmid_type_genes": merlin_magic.abricate_abaum_genes,
                "abricate_abaum_version": merlin_magic.abricate_abaum_version,
                "abricate_database": abricate.abricate_database,
                "abricate_docker": abricate.abricate_docker,
                "abricate_genes": abricate.abricate_genes,
                "abricate_results_tsv": abricate.abricate_results,
                "abricate_version": abricate.abricate_version,
                "abricate_vibrio_biotype": merlin_magic.abricate_vibrio_biotype,
                "abricate_vibrio_ctxA": merlin_magic.abricate_vibrio_ctxA,
                "abricate_vibrio_database": merlin_magic.abricate_vibrio_database,
                "abricate_vibrio_detailed_tsv": merlin_magic.abricate_vibrio_detailed_tsv,
                "abricate_vibrio_docker": merlin_magic.abricate_vibrio_docker,
                "abricate_vibrio_ompW": merlin_magic.abricate_vibrio_ompW,
                "abricate_vibrio_serogroup": merlin_magic.abricate_vibrio_serogroup,
                "abricate_vibrio_toxR": merlin_magic.abricate_vibrio_toxR,
                "abricate_vibrio_version": merlin_magic.abricate_vibrio_version,
                "agrvate_agr_canonical": merlin_magic.agrvate_agr_canonical,
                "agrvate_agr_group": merlin_magic.agrvate_agr_group,
                "agrvate_agr_match_score": merlin_magic.agrvate_agr_match_score,
                "agrvate_agr_multiple": merlin_magic.agrvate_agr_multiple,
                "agrvate_agr_num_frameshifts": merlin_magic.agrvate_agr_num_frameshifts,
                "agrvate_docker": merlin_magic.agrvate_docker,
                "agrvate_results": merlin_magic.agrvate_results,
                "agrvate_summary": merlin_magic.agrvate_summary,
                "agrvate_version": merlin_magic.agrvate_version,    
                "amr_search_csv": merlin_magic.amr_results_csv,
                "amr_search_docker": merlin_magic.amr_search_docker,
                "amr_search_results": merlin_magic.amr_search_results,
                "amr_search_results_pdf": merlin_magic.amr_results_pdf,
                "amr_search_version": merlin_magic.amr_search_version,
                "amrfinderplus_all_report": amrfinderplus_task.amrfinderplus_all_report,
                "amrfinderplus_amr_betalactam_betalactam_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_betalactam_genes,
                "amrfinderplus_amr_betalactam_carbapenem_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_carbapenem_genes,
                "amrfinderplus_amr_betalactam_cephalosporin_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_cephalosporin_genes,
                "amrfinderplus_amr_betalactam_cephalothin_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_cephalothin_genes,
                "amrfinderplus_amr_betalactam_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_genes,
                "amrfinderplus_amr_betalactam_methicillin_genes": amrfinderplus_task.amrfinderplus_amr_betalactam_methicillin_genes,
                "amrfinderplus_amr_classes": amrfinderplus_task.amrfinderplus_amr_classes,
                "amrfinderplus_amr_core_genes": amrfinderplus_task.amrfinderplus_amr_core_genes,
                "amrfinderplus_amr_plus_genes": amrfinderplus_task.amrfinderplus_amr_plus_genes,
                "amrfinderplus_amr_report": amrfinderplus_task.amrfinderplus_amr_report,
                "amrfinderplus_amr_subclasses": amrfinderplus_task.amrfinderplus_amr_subclasses,
                "amrfinderplus_db_version": amrfinderplus_task.amrfinderplus_db_version,
                "amrfinderplus_stress_genes": amrfinderplus_task.amrfinderplus_stress_genes,
                "amrfinderplus_stress_report": amrfinderplus_task.amrfinderplus_stress_report,
                "amrfinderplus_version": amrfinderplus_task.amrfinderplus_version,
                "amrfinderplus_virulence_genes": amrfinderplus_task.amrfinderplus_virulence_genes,
                "amrfinderplus_virulence_report": amrfinderplus_task.amrfinderplus_virulence_report,
                "ani_highest_percent_bases_aligned": ani.ani_highest_percent_bases_aligned,
                "ani_highest_percent": ani.ani_highest_percent,
                "ani_mummer_docker": ani.ani_docker,
                "ani_mummer_version": ani.ani_mummer_version,
                "ani_output_tsv": ani.ani_output_tsv,
                "ani_top_species_match": ani.ani_top_species_match,
                "assembly_fasta": digger_denovo.assembly_fasta,
                "assembly_length": quast.genome_length,
                "bakta_gbff": bakta.bakta_gbff,
                "bakta_gff3": bakta.bakta_gff3,
                "bakta_plot": bakta.bakta_plot,
                "bakta_summary": bakta.bakta_txt,
                "bakta_tsv": bakta.bakta_tsv,
                "bakta_version": bakta.bakta_version,
                "bbduk_docker": read_QC_trim.bbduk_docker,
                "busco_database": busco.busco_database,
                "busco_docker": busco.busco_docker,
                "busco_report": busco.busco_report,
                "busco_results": busco.busco_results,
                "busco_version": busco.busco_version,
                "cg_pipeline_docker": cg_pipeline_raw.cg_pipeline_docker,
                "cg_pipeline_report_clean": cg_pipeline_clean.cg_pipeline_report,
                "cg_pipeline_report_raw": cg_pipeline_raw.cg_pipeline_report,
                "city": city,
                "collection_date": collection_date,
                "contigs_gfa": digger_denovo.contigs_gfa,
                "county": county,
                "ectyper_database_version": merlin_magic.ectyper_database_version,
                "ectyper_docker": merlin_magic.ectyper_docker,
                "ectyper_pathodb_version": merlin_magic.ectyper_pathodb_version,
                "ectyper_pathotype": merlin_magic.ectyper_pathotype,
                "ectyper_pathotype_count": merlin_magic.ectyper_pathotype_count,
                "ectyper_pathotype_genes": merlin_magic.ectyper_pathotype_genes,
                "ectyper_qc_result": merlin_magic.ectyper_qc_result,
                "ectyper_predicted_serotype": merlin_magic.ectyper_predicted_serotype,
                "ectyper_results": merlin_magic.ectyper_results,
                "ectyper_stx_subtypes": merlin_magic.ectyper_stx_subtypes,
                "ectyper_version": merlin_magic.ectyper_version,
                "ectyper_warnings": merlin_magic.ectyper_warnings,
                "emmtyper_docker": merlin_magic.emmtyper_docker,
                "emmtyper_emm_type": merlin_magic.emmtyper_emm_type,
                "emmtyper_results_tsv": merlin_magic.emmtyper_results_tsv,
                "emmtyper_version": merlin_magic.emmtyper_version,
                "est_coverage_clean": cg_pipeline_clean.est_coverage,
                "est_coverage_raw": cg_pipeline_raw.est_coverage,
                "fastp_html_report": read_QC_trim.fastp_html_report,
                "fastp_json_report": read_QC_trim.fastp_json_report,
                "fastp_docker": read_QC_trim.fastp_docker,
                "fastp_version": read_QC_trim.fastp_version,
                "fastq_scan_clean1_json": read_QC_trim.fastq_scan_clean1_json,
                "fastq_scan_num_reads_clean1": read_QC_trim.fastq_scan_clean1,
                "fastq_scan_num_reads_raw1": read_QC_trim.fastq_scan_raw1,
                "fastq_scan_raw1_json": read_QC_trim.fastq_scan_raw1_json,
                "fastq_scan_version": read_QC_trim.fastq_scan_version,
                "fastqc_clean1_html": read_QC_trim.fastqc_clean1_html,
                "fastqc_docker": read_QC_trim.fastqc_docker,
                "fastqc_num_reads_clean1": read_QC_trim.fastqc_clean1,
                "fastqc_num_reads_raw1": read_QC_trim.fastqc_raw1,
                "fastqc_raw1_html": read_QC_trim.fastqc_raw1_html,
                "fastqc_version": read_QC_trim.fastqc_version,
                "gambit_closest_genomes": gambit.gambit_closest_genomes_file,
                "gambit_db_version": gambit.gambit_db_version,
                "gambit_docker": gambit.gambit_docker,
                "gambit_predicted_taxon_rank": gambit.gambit_predicted_taxon_rank,
                "gambit_predicted_taxon": gambit.gambit_predicted_taxon,
                "gambit_report": gambit.gambit_report_file,
                "gambit_version": gambit.gambit_version,
                "genotyphi_final_genotype": merlin_magic.genotyphi_final_genotype,
                "genotyphi_genotype_confidence": merlin_magic.genotyphi_genotype_confidence,
                "genotyphi_mykrobe_json": merlin_magic.genotyphi_mykrobe_json,
                "genotyphi_report_tsv": merlin_magic.genotyphi_report_tsv,
                "genotyphi_species": merlin_magic.genotyphi_species,
                "genotyphi_st_probes_percent_coverage": merlin_magic.genotyphi_st_probes_percent_coverage,
                "genotyphi_version": merlin_magic.genotyphi_version,
                "hicap_docker": merlin_magic.hicap_docker,
                "hicap_genes": merlin_magic.hicap_genes,
                "hicap_results_tsv": merlin_magic.hicap_results_tsv,
                "hicap_serotype": merlin_magic.hicap_serotype,
                "hicap_version": merlin_magic.hicap_version,
                "kaptive_k_locus": merlin_magic.kaptive_k_match,
                "kaptive_k_type": merlin_magic.kaptive_k_type,
                "kaptive_kl_confidence": merlin_magic.kaptive_k_confidence,
                "kaptive_oc_locus": merlin_magic.kaptive_oc_match,
                "kaptive_ocl_confidence": merlin_magic.kaptive_oc_confidence,
                "kaptive_output_file_k": merlin_magic.kaptive_output_file_k,
                "kaptive_output_file_oc": merlin_magic.kaptive_output_file_oc,
                "kaptive_version": merlin_magic.kaptive_version,
                "kleborate_docker": merlin_magic.kleborate_docker,
                "kleborate_genomic_resistance_mutations": merlin_magic.kleborate_genomic_resistance_mutations,
                "kleborate_key_resistance_genes": merlin_magic.kleborate_key_resistance_genes,
                "kleborate_klocus_confidence": merlin_magic.kleborate_klocus_confidence,
                "kleborate_klocus": merlin_magic.kleborate_klocus,
                "kleborate_ktype": merlin_magic.kleborate_ktype,
                "kleborate_mlst_sequence_type": merlin_magic.kleborate_mlst_sequence_type,
                "kleborate_olocus_confidence": merlin_magic.kleborate_olocus_confidence,
                "kleborate_olocus": merlin_magic.kleborate_olocus,
                "kleborate_otype": merlin_magic.kleborate_otype,
                "kleborate_output_file": merlin_magic.kleborate_output_file,
                "kleborate_resistance_score": merlin_magic.kleborate_resistance_score,
                "kleborate_version": merlin_magic.kleborate_version,
                "kleborate_virulence_score": merlin_magic.kleborate_virulence_score,
                "kmerfinder_database": kmerfinder.kmerfinder_database,
                "kmerfinder_docker": kmerfinder.kmerfinder_docker,
                "kmerfinder_query_coverage": kmerfinder.kmerfinder_query_coverage,
                "kmerfinder_results_tsv": kmerfinder.kmerfinder_results_tsv,
                "kmerfinder_template_coverage": kmerfinder.kmerfinder_template_coverage,
                "kmerfinder_top_hit": kmerfinder.kmerfinder_top_hit,
                "kraken2_database": read_QC_trim.kraken_database,
                "kraken2_docker": read_QC_trim.kraken_docker,
                "kraken2_report": read_QC_trim.kraken_report,
                "kraken2_version": read_QC_trim.kraken_version,
                "legsta_predicted_sbt": merlin_magic.legsta_predicted_sbt,
                "legsta_results": merlin_magic.legsta_results,
                "legsta_version": merlin_magic.legsta_version,
                "lissero_results": merlin_magic.lissero_results,
                "lissero_serotype": merlin_magic.lissero_serotype,
                "lissero_version": merlin_magic.lissero_version,
                "meningotype_BAST": merlin_magic.meningotype_BAST,
                "meningotype_FetA": merlin_magic.meningotype_FetA,
                "meningotype_fHbp": merlin_magic.meningotype_fHbp,
                "meningotype_NadA": merlin_magic.meningotype_NadA,
                "meningotype_NHBA": merlin_magic.meningotype_NHBA,
                "meningotype_PorA": merlin_magic.meningotype_PorA,
                "meningotype_PorB": merlin_magic.meningotype_PorB,
                "meningotype_serogroup": merlin_magic.meningotype_serogroup,
                "meningotype_tsv": merlin_magic.meningotype_tsv,
                "meningotype_version": merlin_magic.meningotype_version,
                "midas_docker": read_QC_trim.midas_docker,
                "midas_primary_genus": read_QC_trim.midas_primary_genus,
                "midas_report": read_QC_trim.midas_report,
                "midas_secondary_genus_abundance": read_QC_trim.midas_secondary_genus_abundance,
                "midas_secondary_genus_coverage": read_QC_trim.midas_secondary_genus_coverage,
                "midas_secondary_genus": read_QC_trim.midas_secondary_genus,
                "n50_value": quast.n50_value,
                "ngmaster_ngmast_porB_allele": merlin_magic.ngmaster_ngmast_porB_allele,
                "ngmaster_ngmast_sequence_type": merlin_magic.ngmaster_ngmast_sequence_type,
                "ngmaster_ngmast_tbpB_allele": merlin_magic.ngmaster_ngmast_tbpB_allele,
                "ngmaster_ngstar_23S_allele": merlin_magic.ngmaster_ngstar_23S_allele,
                "ngmaster_ngstar_gyrA_allele": merlin_magic.ngmaster_ngstar_gyrA_allele,
                "ngmaster_ngstar_mtrR_allele": merlin_magic.ngmaster_ngstar_mtrR_allele,
                "ngmaster_ngstar_parC_allele": merlin_magic.ngmaster_ngstar_parC_allele,
                "ngmaster_ngstar_penA_allele": merlin_magic.ngmaster_ngstar_penA_allele,
                "ngmaster_ngstar_ponA_allele": merlin_magic.ngmaster_ngstar_ponA_allele,
                "ngmaster_ngstar_porB_allele": merlin_magic.ngmaster_ngstar_porB_allele,
                "ngmaster_ngstar_sequence_type": merlin_magic.ngmaster_ngstar_sequence_type,
                "ngmaster_tsv": merlin_magic.ngmaster_tsv,
                "ngmaster_version": merlin_magic.ngmaster_version,
                "number_contigs": quast.number_contigs,
                "originating_lab": originating_lab,
                "pasty_all_serogroups": merlin_magic.pasty_all_serogroups,
                "pasty_blast_hits": merlin_magic.pasty_blast_hits,
                "pasty_comment": merlin_magic.pasty_comment,
                "pasty_docker": merlin_magic.pasty_docker,
                "pasty_serogroup_coverage": merlin_magic.pasty_serogroup_coverage,
                "pasty_serogroup_fragments": merlin_magic.pasty_serogroup_fragments,
                "pasty_serogroup": merlin_magic.pasty_serogroup,
                "pasty_summary_tsv": merlin_magic.pasty_summary_tsv,
                "pasty_version": merlin_magic.pasty_version,
                "pbptyper_docker": merlin_magic.pbptyper_docker,
                "pbptyper_pbptype_predicted_tsv": merlin_magic.pbptyper_pbptype_predicted_tsv,
                "pbptyper_predicted_1A_2B_2X": merlin_magic.pbptyper_predicted_1A_2B_2X,
                "pbptyper_version": merlin_magic.pbptyper_version,
                "plasmidfinder_db_version": plasmidfinder.plasmidfinder_db_version,
                "plasmidfinder_docker": plasmidfinder.plasmidfinder_docker,
                "plasmidfinder_plasmids": plasmidfinder.plasmidfinder_plasmids,
                "plasmidfinder_results": plasmidfinder.plasmidfinder_results,
                "plasmidfinder_seqs": plasmidfinder.plasmidfinder_seqs,
                "poppunk_docker": merlin_magic.poppunk_docker,
                "poppunk_gps_cluster": merlin_magic.poppunk_gps_cluster,
                "poppunk_GPS_db_version": merlin_magic.poppunk_GPS_db_version,
                "poppunk_gps_external_cluster_csv": merlin_magic.poppunk_gps_external_cluster_csv,
                "poppunk_version": merlin_magic.poppunk_version,
                "prokka_gbk": prokka.prokka_gbk,
                "prokka_gff": prokka.prokka_gff,
                "prokka_sqn": prokka.prokka_sqn,
                "qc_check": qc_check_task.qc_check,
                "qc_standard": qc_check_task.qc_standard,
                "quast_gc_percent": quast.gc_percent,
                "quast_report": quast.quast_report,
                "quast_version": quast.version,
                "r1_mean_q_clean": cg_pipeline_clean.r1_mean_q,
                "r1_mean_q_raw": cg_pipeline_raw.r1_mean_q,
                "r1_mean_readlength_clean": cg_pipeline_clean.r1_mean_readlength,
                "r1_mean_readlength_raw": cg_pipeline_raw.r1_mean_readlength,
                "read_screen_clean_tsv": clean_check_reads.read_screen_tsv,
                "read_screen_clean": clean_check_reads.read_screen,
                "read_screen_raw_tsv": raw_check_reads.read_screen_tsv,
                "read_screen_raw": raw_check_reads.read_screen,
                "read1_clean": read_QC_trim.read1_clean,
                "read1_concatenated": concatenate_illumina_lanes.read1_concatenated,
                "read1": select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
                "resfinder_db_version": resfinder_task.resfinder_db_version,
                "resfinder_docker": resfinder_task.resfinder_docker,
                "resfinder_pheno_table_species": resfinder_task.resfinder_pheno_table_species,
                "resfinder_pheno_table": resfinder_task.resfinder_pheno_table,
                "resfinder_pointfinder_pheno_table": resfinder_task.pointfinder_pheno_table,
                "resfinder_pointfinder_results": resfinder_task.pointfinder_results,
                "resfinder_predicted_pheno_resistance": resfinder_task.resfinder_predicted_pheno_resistance,
                "resfinder_predicted_resistance_Amp": resfinder_task.resfinder_predicted_resistance_Amp,
                "resfinder_predicted_resistance_Axo": resfinder_task.resfinder_predicted_resistance_Axo,
                "resfinder_predicted_resistance_Azm": resfinder_task.resfinder_predicted_resistance_Azm,
                "resfinder_predicted_resistance_Cip": resfinder_task.resfinder_predicted_resistance_Cip,
                "resfinder_predicted_resistance_Smx": resfinder_task.resfinder_predicted_resistance_Smx,
                "resfinder_predicted_resistance_Tmp": resfinder_task.resfinder_predicted_resistance_Tmp,
                "resfinder_predicted_resistance_quinolone": resfinder_task.resfinder_predicted_resistance_quinolone,
                "resfinder_predicted_resistance_quinolone_mechanisms": resfinder_task.resfinder_predicted_resistance_quinolone_mechanisms,
                "resfinder_predicted_xdr_shigella": resfinder_task.resfinder_predicted_xdr_shigella,
                "resfinder_results": resfinder_task.resfinder_results_tab,
                "resfinder_seqs": resfinder_task.resfinder_hit_in_genome_seq,
                "resfinder_version": resfinder_task.resfinder_version,
                "run_id": run_id,
                "seq_platform": seq_method,
                "seqsero2s_note": merlin_magic.seqsero2s_note,
                "seqsero2s_predicted_antigenic_profile": merlin_magic.seqsero2s_predicted_antigenic_profile,
                "seqsero2s_predicted_contamination": merlin_magic.seqsero2s_predicted_contamination,
                "seqsero2s_predicted_serotype": merlin_magic.seqsero2s_predicted_serotype,
                "seqsero2s_report": merlin_magic.seqsero2s_report,
                "seqsero2s_version": merlin_magic.seqsero2s_version,
                "serotypefinder_docker": merlin_magic.serotypefinder_docker,
                "serotypefinder_report": merlin_magic.serotypefinder_report,
                "serotypefinder_serotype": merlin_magic.serotypefinder_serotype,
                "shigatyper_docker": merlin_magic.shigatyper_docker,
                "shigatyper_hits_tsv": merlin_magic.shigatyper_hits_tsv,
                "shigatyper_ipaB_presence_absence": merlin_magic.shigatyper_ipaB_presence_absence,
                "shigatyper_notes": merlin_magic.shigatyper_notes,
                "shigatyper_predicted_serotype": merlin_magic.shigatyper_predicted_serotype,
                "shigatyper_summary_tsv": merlin_magic.shigatyper_summary_tsv,
                "shigatyper_version": merlin_magic.shigatyper_version,
                "shigeifinder_cluster_reads": merlin_magic.shigeifinder_cluster_reads,
                "shigeifinder_cluster": merlin_magic.shigeifinder_cluster,
                "shigeifinder_docker_reads": merlin_magic.shigeifinder_docker_reads,
                "shigeifinder_docker": merlin_magic.shigeifinder_docker,
                "shigeifinder_H_antigen_reads": merlin_magic.shigeifinder_H_antigen_reads,
                "shigeifinder_H_antigen": merlin_magic.shigeifinder_H_antigen,
                "shigeifinder_ipaH_presence_absence_reads": merlin_magic.shigeifinder_ipaH_presence_absence_reads,
                "shigeifinder_ipaH_presence_absence": merlin_magic.shigeifinder_ipaH_presence_absence,
                "shigeifinder_notes_reads": merlin_magic.shigeifinder_notes_reads,
                "shigeifinder_notes": merlin_magic.shigeifinder_notes,
                "shigeifinder_num_virulence_plasmid_genes_reads": merlin_magic.shigeifinder_num_virulence_plasmid_genes_reads,
                "shigeifinder_num_virulence_plasmid_genes": merlin_magic.shigeifinder_num_virulence_plasmid_genes,
                "shigeifinder_O_antigen_reads": merlin_magic.shigeifinder_O_antigen_reads,
                "shigeifinder_O_antigen": merlin_magic.shigeifinder_O_antigen,
                "shigeifinder_report_reads": merlin_magic.shigeifinder_report_reads,
                "shigeifinder_report": merlin_magic.shigeifinder_report,
                "shigeifinder_serotype_reads": merlin_magic.shigeifinder_serotype_reads,
                "shigeifinder_serotype": merlin_magic.shigeifinder_serotype,
                "shigeifinder_version_reads": merlin_magic.shigeifinder_version_reads,
                "shigeifinder_version": merlin_magic.shigeifinder_version,
                "assembler_version": digger_denovo.assembler_version,
                "sistr_allele_fasta": merlin_magic.sistr_allele_fasta,
                "sistr_allele_json": merlin_magic.sistr_allele_json,
                "sistr_antigenic_formula": merlin_magic.sistr_antigenic_formula,
                "sistr_cgmlst": merlin_magic.sistr_cgmlst,
                "sistr_h1_antigens": merlin_magic.sistr_h1_antigens,
                "sistr_h2_antigens": merlin_magic.sistr_h2_antigens,
                "sistr_o_antigens": merlin_magic.sistr_o_antigens,
                "sistr_predicted_serotype": merlin_magic.sistr_predicted_serotype,
                "sistr_results": merlin_magic.sistr_results,                
                "sistr_serogroup": merlin_magic.sistr_serogroup,
                "sistr_serotype_cgmlst": merlin_magic.sistr_serotype_cgmlst,
                "sistr_version": merlin_magic.sistr_version,
                "sonneityping_final_genotype": merlin_magic.sonneityping_final_genotype,
                "sonneityping_final_report_tsv": merlin_magic.sonneityping_final_report_tsv,
                "sonneityping_genotype_confidence": merlin_magic.sonneityping_genotype_confidence,
                "sonneityping_genotype_name": merlin_magic.sonneityping_genotype_name,
                "sonneityping_mykrobe_docker": merlin_magic.sonneityping_mykrobe_docker,
                "sonneityping_mykrobe_report_csv": merlin_magic.sonneityping_mykrobe_report_csv,
                "sonneityping_mykrobe_report_json": merlin_magic.sonneityping_mykrobe_report_json,
                "sonneityping_mykrobe_version": merlin_magic.sonneityping_mykrobe_version,
                "sonneityping_species": merlin_magic.sonneityping_species,
                "spatyper_docker": merlin_magic.spatyper_docker,
                "spatyper_repeats": merlin_magic.spatyper_repeats,
                "spatyper_tsv": merlin_magic.spatyper_tsv,
                "spatyper_type": merlin_magic.spatyper_type,
                "spatyper_version": merlin_magic.spatyper_version,
                "srst2_vibrio_biotype": merlin_magic.srst2_vibrio_biotype,
                "srst2_vibrio_ctxA": merlin_magic.srst2_vibrio_ctxA,
                "srst2_vibrio_database": merlin_magic.srst2_vibrio_database,
                "srst2_vibrio_detailed_tsv": merlin_magic.srst2_vibrio_detailed_tsv,
                "srst2_vibrio_docker": merlin_magic.srst2_vibrio_docker,
                "srst2_vibrio_ompW": merlin_magic.srst2_vibrio_ompW,
                "srst2_vibrio_serogroup": merlin_magic.srst2_vibrio_serogroup,
                "srst2_vibrio_toxR": merlin_magic.srst2_vibrio_toxR,
                "srst2_vibrio_version": merlin_magic.srst2_vibrio_version,
                "staphopiasccmec_docker": merlin_magic.staphopiasccmec_docker,
                "staphopiasccmec_hamming_distance_tsv": merlin_magic.staphopiasccmec_hamming_distance_tsv,
                "staphopiasccmec_results_tsv": merlin_magic.staphopiasccmec_results_tsv,
                "staphopiasccmec_types_and_mecA_presence": merlin_magic.staphopiasccmec_types_and_mecA_presence,
                "staphopiasccmec_version": merlin_magic.staphopiasccmec_version,
                "stxtyper_all_hits": merlin_magic.stxtyper_all_hits,
                "stxtyper_ambiguous_hits": merlin_magic.stxtyper_ambiguous_hits,
                "stxtyper_complete_operons": merlin_magic.stxtyper_complete_operon_hits,
                "stxtyper_docker": merlin_magic.stxtyper_docker,
                "stxtyper_extended_operons": merlin_magic.stxtyper_extended_operons,
                "stxtyper_novel_hits": merlin_magic.stxtyper_novel_hits,
                "stxtyper_num_hits": merlin_magic.stxtyper_num_hits,
                "stxtyper_partial_hits": merlin_magic.stxtyper_partial_hits,
                "stxtyper_report": merlin_magic.stxtyper_report,
                "stxtyper_stx_frameshifts_or_internal_stop_hits": merlin_magic.stxtyper_stx_frameshifts_or_internal_stop_hits,
                "stxtyper_version": merlin_magic.stxtyper_version,
                "tbp_parser_average_genome_depth": merlin_magic.tbp_parser_average_genome_depth,
                "tbp_parser_coverage_report": merlin_magic.tbp_parser_coverage_report,
                "tbp_parser_genome_percent_coverage": merlin_magic.tbp_parser_genome_percent_coverage,
                "tbp_parser_laboratorian_report_csv": merlin_magic.tbp_parser_laboratorian_report_csv,
                "tbp_parser_lims_report_csv": merlin_magic.tbp_parser_lims_report_csv,
                "tbp_parser_looker_report_csv": merlin_magic.tbp_parser_looker_report_csv,
                "tbprofiler_dr_type": merlin_magic.tbprofiler_dr_type,
                "tbprofiler_main_lineage": merlin_magic.tbprofiler_main_lineage,
                "tbprofiler_output_bai": merlin_magic.tbprofiler_output_bai,
                "tbprofiler_output_bam": merlin_magic.tbprofiler_output_bam,
                "tbprofiler_output_file": merlin_magic.tbprofiler_output_file,
                "tbprofiler_output_vcf": merlin_magic.tbprofiler_output_vcf,
                "tbprofiler_resistance_genes": merlin_magic.tbprofiler_resistance_genes,
                "tbprofiler_sub_lineage": merlin_magic.tbprofiler_sub_lineage,
                "tbprofiler_version": merlin_magic.tbprofiler_version,
                "theiaprok_illumina_se_analysis_date": version_capture.date,
                "theiaprok_illumina_se_version": version_capture.phb_version,
                "trimmomatic_docker": read_QC_trim.trimmomatic_docker,
                "trimmomatic_version": read_QC_trim.trimmomatic_version,
                "ts_mlst_allelic_profile": ts_mlst.ts_mlst_allelic_profile,
                "ts_mlst_docker": ts_mlst.ts_mlst_docker,
                "ts_mlst_novel_alleles": ts_mlst.ts_mlst_novel_alleles,
                "ts_mlst_predicted_st": ts_mlst.ts_mlst_predicted_st,
                "ts_mlst_pubmlst_scheme": ts_mlst.ts_mlst_pubmlst_scheme,
                "ts_mlst_predicted_secondary_st": ts_mlst.ts_mlst_predicted_secondary_st,
                "ts_mlst_pubmlst_secondary_scheme": ts_mlst.ts_mlst_pubmlst_secondary_scheme,
                "ts_mlst_secondary_allelic_profile": ts_mlst.ts_mlst_secondary_allelic_profile,
                "ts_mlst_secondary_novel_alleles": ts_mlst.ts_mlst_secondary_novel_alleles,
                "ts_mlst_results": ts_mlst.ts_mlst_results,
                "ts_mlst_version": ts_mlst.ts_mlst_version,
                "virulencefinder_docker": merlin_magic.virulencefinder_docker,
                "virulencefinder_hits": merlin_magic.virulencefinder_hits,
                "virulencefinder_report_tsv": merlin_magic.virulencefinder_report_tsv,
                "zip": zip,
                "arln_assembly_ratio": arln_stats.assembly_ratio,
                "arln_assembly_zscore": arln_stats.assembly_zscore,
                "arln_r1_q30_clean": arln_stats.read1_clean_q30,
                "arln_r1_q30_raw": arln_stats.read1_raw_q30,
                "arln_stats_docker_version": arln_stats.docker_version,
                "arln_taxon_assembly_ratio_stdev": arln_stats.taxon_assembly_ratio_stdev,
                "arln_taxon_gc_mean": arln_stats.taxon_gc_mean,
                "arln_taxon_gc_percent_stdev": arln_stats.taxon_gc_percent_stdev,
                "assembler": digger_denovo.assembler_used,
                "filtered_contigs_metrics": digger_denovo.filtered_contigs_metrics,
                "gamma_docker": gamma.gamma_docker,
                "gamma_fasta": gamma.gamma_fasta,
                "gamma_gff": gamma.gamma_gff,
                "gamma_results": gamma.gamma_results,
                "gamma_version": gamma.gamma_version
            }
          }
        }
        if (call_arln_stats) {
          call arln_stats.arln_stats {
            input:
              samplename = samplename,
              taxon = select_first([gambit.gambit_predicted_taxon, expected_taxon]),
              workflow_type = "se",
              genome_length = quast.genome_length,
              gc_percent = quast.gc_percent,
              read1_raw = select_first([concatenate_illumina_lanes.read1_concatenated, read1]),
              read1_clean = read_QC_trim.read1_clean
          }
        }
      }
    }
  }
  output {
    # Version Captures
    String theiaprok_illumina_se_version = version_capture.phb_version
    String theiaprok_illumina_se_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Concatenated Illumina Reads
    File? read1_concatenated = concatenate_illumina_lanes.read1_concatenated
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC - fastq_scan outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # Read QC - bbduk outputs
    File? read1_clean = read_QC_trim.read1_clean
    String? bbduk_docker = read_QC_trim.bbduk_docker
    # Read QC - cg pipeline outputs
    Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
    Float? r1_mean_q_clean = cg_pipeline_clean.r1_mean_q
    Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
    Float? r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength
    # Read QC - midas outputs
    String? midas_docker = read_QC_trim.midas_docker
    File? midas_report = read_QC_trim.midas_report
    String? midas_primary_genus = read_QC_trim.midas_primary_genus
    String? midas_secondary_genus = read_QC_trim.midas_secondary_genus
    Float? midas_secondary_genus_abundance = read_QC_trim.midas_secondary_genus_abundance
    Float? midas_secondary_genus_coverage = read_QC_trim.midas_secondary_genus_coverage
    # Read QC - kraken outputs
    String? kraken2_version = read_QC_trim.kraken_version
    String? kraken2_report = read_QC_trim.kraken_report
    String? kraken2_docker = read_QC_trim.kraken_docker
    String? kraken2_database = read_QC_trim.kraken_database
    #Assembly - digger_denovo outputs
    File? assembly_fasta = digger_denovo.assembly_fasta
    File? contigs_gfa = digger_denovo.contigs_gfa
    File? filtered_contigs_metrics = digger_denovo.filtered_contigs_metrics
    String? assembler = digger_denovo.assembler_used
    String? assembler_version = digger_denovo.assembler_version
    # Assembly QC - quast outputs
    File? quast_report = quast.quast_report
    String? quast_version = quast.version
    Int? assembly_length = quast.genome_length
    Int? number_contigs = quast.number_contigs
    Int? n50_value = quast.n50_value
    Float? quast_gc_percent = quast.gc_percent
    # Assembly QC - cg pipeline outputs
    File? cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report
    String? cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker
    Float? est_coverage_raw = cg_pipeline_raw.est_coverage
    File? cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report
    Float? est_coverage_clean = cg_pipeline_clean.est_coverage
    # Assembly QC - busco outputs
    String? busco_version = busco.busco_version
    String? busco_docker = busco.busco_docker
    String? busco_database = busco.busco_database
    String? busco_results = busco.busco_results
    File? busco_report = busco.busco_report
    # Taxon ID - gambit outputs
    File? gambit_report = gambit.gambit_report_file
    File? gambit_closest_genomes = gambit.gambit_closest_genomes_file
    String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String? gambit_version = gambit.gambit_version
    String? gambit_db_version = gambit.gambit_db_version
    String? gambit_docker = gambit.gambit_docker
    # ani-mummer outputs
    Float? ani_highest_percent = ani.ani_highest_percent
    Float? ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned
    File? ani_output_tsv = ani.ani_output_tsv
    String? ani_top_species_match = ani.ani_top_species_match
    String? ani_mummer_version = ani.ani_mummer_version
    String? ani_mummer_docker = ani.ani_docker
    # kmerfinder outputs
    String? kmerfinder_docker = kmerfinder.kmerfinder_docker
    File? kmerfinder_results_tsv = kmerfinder.kmerfinder_results_tsv
    String? kmerfinder_top_hit = kmerfinder.kmerfinder_top_hit
    String? kmerfinder_query_coverage = kmerfinder.kmerfinder_query_coverage
    String? kmerfinder_template_coverage = kmerfinder.kmerfinder_template_coverage
    String? kmerfinder_database = kmerfinder.kmerfinder_database
    # NCBI-AMRFinderPlus Outputs
    File? amrfinderplus_all_report = amrfinderplus_task.amrfinderplus_all_report
    File? amrfinderplus_amr_report = amrfinderplus_task.amrfinderplus_amr_report
    File? amrfinderplus_stress_report = amrfinderplus_task.amrfinderplus_stress_report
    File? amrfinderplus_virulence_report = amrfinderplus_task.amrfinderplus_virulence_report
    String? amrfinderplus_amr_core_genes = amrfinderplus_task.amrfinderplus_amr_core_genes
    String? amrfinderplus_amr_plus_genes = amrfinderplus_task.amrfinderplus_amr_plus_genes
    String? amrfinderplus_stress_genes = amrfinderplus_task.amrfinderplus_stress_genes
    String? amrfinderplus_virulence_genes = amrfinderplus_task.amrfinderplus_virulence_genes
    String? amrfinderplus_amr_classes = amrfinderplus_task.amrfinderplus_amr_classes
    String? amrfinderplus_amr_subclasses = amrfinderplus_task.amrfinderplus_amr_subclasses
    String? amrfinderplus_version = amrfinderplus_task.amrfinderplus_version
    String? amrfinderplus_db_version = amrfinderplus_task.amrfinderplus_db_version
    # NCBI-AMRFinderPlus Outputs for BETA-LACTAM genes
    String? amrfinderplus_amr_betalactam_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_genes
    String? amrfinderplus_amr_betalactam_betalactam_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_betalactam_genes
    String? amrfinderplus_amr_betalactam_carbapenem_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_carbapenem_genes
    String? amrfinderplus_amr_betalactam_cephalosporin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_cephalosporin_genes
    String? amrfinderplus_amr_betalactam_cephalothin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_cephalothin_genes
    String? amrfinderplus_amr_betalactam_methicillin_genes = amrfinderplus_task.amrfinderplus_amr_betalactam_methicillin_genes
    # GAMMA Outputs
    File? gamma_results = gamma.gamma_results
    File? gamma_gff = gamma.gamma_gff
    File? gamma_fasta = gamma.gamma_fasta
    String? gamma_version = gamma.gamma_version
    String? gamma_docker = gamma.gamma_docker
    # AMR_Search
    File? amr_search_results = merlin_magic.amr_search_results
    File? amr_search_csv = merlin_magic.amr_results_csv
    File? amr_search_results_pdf = merlin_magic.amr_results_pdf
    String? amr_search_docker = merlin_magic.amr_search_docker
    String? amr_search_version = merlin_magic.amr_search_version
    # Resfinder Outputs
    File? resfinder_pheno_table = resfinder_task.resfinder_pheno_table
    File? resfinder_pheno_table_species = resfinder_task.resfinder_pheno_table_species
    File? resfinder_seqs = resfinder_task.resfinder_hit_in_genome_seq
    File? resfinder_results = resfinder_task.resfinder_results_tab
    File? resfinder_pointfinder_pheno_table = resfinder_task.pointfinder_pheno_table
    File? resfinder_pointfinder_results = resfinder_task.pointfinder_results
    String? resfinder_predicted_pheno_resistance = resfinder_task.resfinder_predicted_pheno_resistance
    String? resfinder_predicted_xdr_shigella = resfinder_task.resfinder_predicted_xdr_shigella
    String? resfinder_predicted_resistance_Amp = resfinder_task.resfinder_predicted_resistance_Amp
    String? resfinder_predicted_resistance_Azm = resfinder_task.resfinder_predicted_resistance_Azm
    String? resfinder_predicted_resistance_Axo = resfinder_task.resfinder_predicted_resistance_Axo
    String? resfinder_predicted_resistance_Cip = resfinder_task.resfinder_predicted_resistance_Cip
    String? resfinder_predicted_resistance_Smx = resfinder_task.resfinder_predicted_resistance_Smx
    String? resfinder_predicted_resistance_Tmp = resfinder_task.resfinder_predicted_resistance_Tmp
    String? resfinder_predicted_resistance_quinolone = resfinder_task.resfinder_predicted_resistance_quinolone
    Int? resfinder_predicted_resistance_quinolone_mechanisms = resfinder_task.resfinder_predicted_resistance_quinolone_mechanisms
    String? resfinder_version = resfinder_task.resfinder_version
    String? resfinder_db_version = resfinder_task.resfinder_db_version
    String? resfinder_docker = resfinder_task.resfinder_docker
    # MLST Typing
    File? ts_mlst_results = ts_mlst.ts_mlst_results
    String? ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st
    String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
    String? ts_mlst_allelic_profile = ts_mlst.ts_mlst_allelic_profile
    File? ts_mlst_novel_alleles = ts_mlst.ts_mlst_novel_alleles
    String? ts_mlst_predicted_secondary_st = ts_mlst.ts_mlst_predicted_secondary_st
    String? ts_mlst_pubmlst_secondary_scheme = ts_mlst.ts_mlst_pubmlst_secondary_scheme
    String? ts_mlst_secondary_allelic_profile = ts_mlst.ts_mlst_secondary_allelic_profile
    File? ts_mlst_secondary_novel_alleles = ts_mlst.ts_mlst_secondary_novel_alleles
    String? ts_mlst_version = ts_mlst.ts_mlst_version
    String? ts_mlst_docker = ts_mlst.ts_mlst_docker
    # Prokka Results
    File? prokka_gff = prokka.prokka_gff
    File? prokka_gbk = prokka.prokka_gbk
    File? prokka_sqn = prokka.prokka_sqn
    # Bakta Results
    File? bakta_gbff = bakta.bakta_gbff
    File? bakta_gff3 = bakta.bakta_gff3
    File? bakta_tsv = bakta.bakta_tsv
    File? bakta_summary = bakta.bakta_txt
    File? bakta_plot = bakta.bakta_plot
    String? bakta_version = bakta.bakta_version
    # Plasmidfinder Results
    String? plasmidfinder_plasmids = plasmidfinder.plasmidfinder_plasmids
    File? plasmidfinder_results = plasmidfinder.plasmidfinder_results
    File? plasmidfinder_seqs = plasmidfinder.plasmidfinder_seqs
    String? plasmidfinder_docker = plasmidfinder.plasmidfinder_docker
    String? plasmidfinder_db_version = plasmidfinder.plasmidfinder_db_version 
    # Abricate Results
    File? abricate_results_tsv = abricate.abricate_results
    String? abricate_genes = abricate.abricate_genes
    String? abricate_database = abricate.abricate_database
    String? abricate_version = abricate.abricate_version
    String? abricate_docker = abricate.abricate_docker  
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # Ecoli Typing
    File? serotypefinder_report = merlin_magic.serotypefinder_report
    String? serotypefinder_docker = merlin_magic.serotypefinder_docker
    String? serotypefinder_serotype = merlin_magic.serotypefinder_serotype
    File? ectyper_results = merlin_magic.ectyper_results
    String? ectyper_version = merlin_magic.ectyper_version
    File? ectyper_warnings = merlin_magic.ectyper_warnings
    String? ectyper_predicted_serotype = merlin_magic.ectyper_predicted_serotype
    String? ectyper_qc_result = merlin_magic.ectyper_qc_result
    String? ectyper_database_version = merlin_magic.ectyper_database_version
    String? ectyper_pathotype = merlin_magic.ectyper_pathotype
    String? ectyper_pathotype_count = merlin_magic.ectyper_pathotype_count
    String? ectyper_pathotype_genes = merlin_magic.ectyper_pathotype_genes
    String? ectyper_pathodb_version = merlin_magic.ectyper_pathodb_version
    String? ectyper_stx_subtypes = merlin_magic.ectyper_stx_subtypes
    String? ectyper_docker = merlin_magic.ectyper_docker
    String? shigatyper_predicted_serotype = merlin_magic.shigatyper_predicted_serotype
    String? shigatyper_ipaB_presence_absence = merlin_magic.shigatyper_ipaB_presence_absence
    String? shigatyper_notes = merlin_magic.shigatyper_notes
    File? shigatyper_hits_tsv = merlin_magic.shigatyper_hits_tsv
    File? shigatyper_summary_tsv = merlin_magic.shigatyper_summary_tsv
    String? shigatyper_version = merlin_magic.shigatyper_version
    String? shigatyper_docker = merlin_magic.shigatyper_docker
    File? shigeifinder_report = merlin_magic.shigeifinder_report
    String? shigeifinder_docker = merlin_magic.shigeifinder_docker
    String? shigeifinder_version = merlin_magic.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence = merlin_magic.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes = merlin_magic.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster = merlin_magic.shigeifinder_cluster
    String? shigeifinder_serotype = merlin_magic.shigeifinder_serotype
    String? shigeifinder_O_antigen = merlin_magic.shigeifinder_O_antigen
    String? shigeifinder_H_antigen = merlin_magic.shigeifinder_H_antigen
    String? shigeifinder_notes = merlin_magic.shigeifinder_notes
    # ShigeiFinder outputs but for task that uses reads instead of assembly as input
    File? shigeifinder_report_reads = merlin_magic.shigeifinder_report
    String? shigeifinder_docker_reads = merlin_magic.shigeifinder_docker
    String? shigeifinder_version_reads = merlin_magic.shigeifinder_version
    String? shigeifinder_ipaH_presence_absence_reads = merlin_magic.shigeifinder_ipaH_presence_absence
    String? shigeifinder_num_virulence_plasmid_genes_reads = merlin_magic.shigeifinder_num_virulence_plasmid_genes
    String? shigeifinder_cluster_reads = merlin_magic.shigeifinder_cluster
    String? shigeifinder_serotype_reads = merlin_magic.shigeifinder_serotype
    String? shigeifinder_O_antigen_reads = merlin_magic.shigeifinder_O_antigen
    String? shigeifinder_H_antigen_reads = merlin_magic.shigeifinder_H_antigen
    String? shigeifinder_notes_reads = merlin_magic.shigeifinder_notes
    # E coli only typing
    File? virulencefinder_report_tsv = merlin_magic.virulencefinder_report_tsv
    String? virulencefinder_docker = merlin_magic.virulencefinder_docker
    String? virulencefinder_hits = merlin_magic.virulencefinder_hits
    # stxtyper 
    File? stxtyper_report = merlin_magic.stxtyper_report
    String? stxtyper_docker = merlin_magic.stxtyper_docker
    String? stxtyper_version = merlin_magic.stxtyper_version
    Int? stxtyper_num_hits = merlin_magic.stxtyper_num_hits
    String? stxtyper_all_hits = merlin_magic.stxtyper_all_hits
    String? stxtyper_complete_operons = merlin_magic.stxtyper_complete_operon_hits
    String? stxtyper_partial_hits = merlin_magic.stxtyper_partial_hits
    String? stxtyper_stx_frameshifts_or_internal_stop_hits =  merlin_magic.stxtyper_stx_frameshifts_or_internal_stop_hits
    String? stxtyper_novel_hits = merlin_magic.stxtyper_novel_hits
    String? stxtyper_extended_operons = merlin_magic.stxtyper_extended_operons
    String? stxtyper_ambiguous_hits = merlin_magic.stxtyper_ambiguous_hits
    # Shigella sonnei Typing
    File? sonneityping_mykrobe_report_csv = merlin_magic.sonneityping_mykrobe_report_csv
    File? sonneityping_mykrobe_report_json = merlin_magic.sonneityping_mykrobe_report_json
    File? sonneityping_final_report_tsv = merlin_magic.sonneityping_final_report_tsv
    String? sonneityping_mykrobe_version = merlin_magic.sonneityping_mykrobe_version
    String? sonneityping_mykrobe_docker = merlin_magic.sonneityping_mykrobe_docker
    String? sonneityping_species = merlin_magic.sonneityping_species
    String? sonneityping_final_genotype = merlin_magic.sonneityping_final_genotype
    String? sonneityping_genotype_confidence = merlin_magic.sonneityping_genotype_confidence
    String? sonneityping_genotype_name = merlin_magic.sonneityping_genotype_name
    # Listeria Typing
    File? lissero_results = merlin_magic.lissero_results
    String? lissero_version = merlin_magic.lissero_version
    String? lissero_serotype = merlin_magic.lissero_serotype
    # Pseudomonas Aeruginosa Typing
    String? pasty_serogroup = merlin_magic.pasty_serogroup
    Float? pasty_serogroup_coverage = merlin_magic.pasty_serogroup_coverage
    Int? pasty_serogroup_fragments = merlin_magic.pasty_serogroup_fragments
    File? pasty_summary_tsv = merlin_magic.pasty_summary_tsv
    File? pasty_blast_hits = merlin_magic.pasty_blast_hits
    File? pasty_all_serogroups = merlin_magic.pasty_all_serogroups
    String? pasty_version = merlin_magic.pasty_version
    String? pasty_docker = merlin_magic.pasty_docker
    String? pasty_comment = merlin_magic.pasty_comment
    # Salmonella Typing
    File? sistr_results = merlin_magic.sistr_results
    File? sistr_allele_json = merlin_magic.sistr_allele_json
    File? sistr_allele_fasta = merlin_magic.sistr_allele_fasta
    File? sistr_cgmlst = merlin_magic.sistr_cgmlst
    String? sistr_version = merlin_magic.sistr_version
    String? sistr_antigenic_formula = merlin_magic.sistr_antigenic_formula
    String? sistr_predicted_serotype = merlin_magic.sistr_predicted_serotype
    String? sistr_serogroup = merlin_magic.sistr_serogroup
    String? sistr_h1_antigens = merlin_magic.sistr_h1_antigens
    String? sistr_h2_antigens = merlin_magic.sistr_h2_antigens
    String? sistr_o_antigens = merlin_magic.sistr_o_antigens
    String? sistr_serotype_cgmlst = merlin_magic.sistr_serotype_cgmlst
    String? seqsero2s_report = merlin_magic.seqsero2s_report
    String? seqsero2s_version = merlin_magic.seqsero2s_version
    String? seqsero2s_note = merlin_magic.seqsero2s_note
    String? seqsero2s_predicted_antigenic_profile = merlin_magic.seqsero2s_predicted_antigenic_profile
    String? seqsero2s_predicted_serotype = merlin_magic.seqsero2s_predicted_serotype
    String? seqsero2s_predicted_contamination = merlin_magic.seqsero2s_predicted_contamination
    # Salmonella serotype Typhi Typing
    File? genotyphi_report_tsv = merlin_magic.genotyphi_report_tsv 
    File? genotyphi_mykrobe_json = merlin_magic.genotyphi_mykrobe_json
    String? genotyphi_version = merlin_magic.genotyphi_version
    String? genotyphi_species = merlin_magic.genotyphi_species
    Float? genotyphi_st_probes_percent_coverage = merlin_magic.genotyphi_st_probes_percent_coverage
    String? genotyphi_final_genotype = merlin_magic.genotyphi_final_genotype
    String? genotyphi_genotype_confidence = merlin_magic.genotyphi_genotype_confidence
    # Klebsiella Typing
    File? kleborate_output_file = merlin_magic.kleborate_output_file
    String? kleborate_version = merlin_magic.kleborate_version
    String? kleborate_docker = merlin_magic.kleborate_docker
    String? kleborate_key_resistance_genes = merlin_magic.kleborate_key_resistance_genes
    String? kleborate_genomic_resistance_mutations = merlin_magic.kleborate_genomic_resistance_mutations
    String? kleborate_mlst_sequence_type = merlin_magic.kleborate_mlst_sequence_type
    String? kleborate_klocus = merlin_magic.kleborate_klocus
    String? kleborate_ktype = merlin_magic.kleborate_ktype
    String? kleborate_olocus = merlin_magic.kleborate_olocus
    String? kleborate_otype = merlin_magic.kleborate_otype
    String? kleborate_klocus_confidence = merlin_magic.kleborate_klocus_confidence
    String? kleborate_olocus_confidence = merlin_magic.kleborate_olocus_confidence
    String? kleborate_virulence_score = merlin_magic.kleborate_virulence_score
    String? kleborate_resistance_score = merlin_magic.kleborate_resistance_score
    # Neisseria gonorrhoeae Typing
    File? ngmaster_tsv = merlin_magic.ngmaster_tsv
    String? ngmaster_version = merlin_magic.ngmaster_version
    String? ngmaster_ngmast_sequence_type = merlin_magic.ngmaster_ngmast_sequence_type
    String? ngmaster_ngmast_porB_allele = merlin_magic.ngmaster_ngmast_porB_allele
    String? ngmaster_ngmast_tbpB_allele = merlin_magic.ngmaster_ngmast_tbpB_allele
    String? ngmaster_ngstar_sequence_type = merlin_magic.ngmaster_ngstar_sequence_type
    String? ngmaster_ngstar_penA_allele = merlin_magic.ngmaster_ngstar_penA_allele
    String? ngmaster_ngstar_mtrR_allele = merlin_magic.ngmaster_ngstar_mtrR_allele
    String? ngmaster_ngstar_porB_allele = merlin_magic.ngmaster_ngstar_porB_allele
    String? ngmaster_ngstar_ponA_allele = merlin_magic.ngmaster_ngstar_ponA_allele
    String? ngmaster_ngstar_gyrA_allele = merlin_magic.ngmaster_ngstar_gyrA_allele
    String? ngmaster_ngstar_parC_allele = merlin_magic.ngmaster_ngstar_parC_allele
    String? ngmaster_ngstar_23S_allele = merlin_magic.ngmaster_ngstar_23S_allele
    # Neisseria meningitidis Typing
    File? meningotype_tsv = merlin_magic.meningotype_tsv
    String? meningotype_version = merlin_magic.meningotype_version
    String? meningotype_serogroup = merlin_magic.meningotype_serogroup
    String? meningotype_PorA = merlin_magic.meningotype_PorA
    String? meningotype_FetA = merlin_magic.meningotype_FetA
    String? meningotype_PorB = merlin_magic.meningotype_PorB
    String? meningotype_fHbp = merlin_magic.meningotype_fHbp
    String? meningotype_NHBA = merlin_magic.meningotype_NHBA
    String? meningotype_NadA = merlin_magic.meningotype_NadA
    String? meningotype_BAST = merlin_magic.meningotype_BAST
    # Acinetobacter Typing
    File? kaptive_output_file_k = merlin_magic.kaptive_output_file_k
    File? kaptive_output_file_oc = merlin_magic.kaptive_output_file_oc
    String? kaptive_version = merlin_magic.kaptive_version
    String? kaptive_k_locus = merlin_magic.kaptive_k_match
    String? kaptive_k_type = merlin_magic.kaptive_k_type
    String? kaptive_kl_confidence = merlin_magic.kaptive_k_confidence
    String? kaptive_oc_locus = merlin_magic.kaptive_oc_match
    String? kaptive_ocl_confidence = merlin_magic.kaptive_oc_confidence
    File? abricate_abaum_plasmid_tsv = merlin_magic.abricate_abaum_results
    String? abricate_abaum_plasmid_type_genes = merlin_magic.abricate_abaum_genes
    String? abricate_abaum_database = merlin_magic.abricate_abaum_database
    String? abricate_abaum_version = merlin_magic.abricate_abaum_version
    String? abricate_abaum_docker = merlin_magic.abricate_abaum_docker
    # Mycobacterium Typing
    File? tbprofiler_output_file = merlin_magic.tbprofiler_output_file
    File? tbprofiler_output_bam = merlin_magic.tbprofiler_output_bam
    File? tbprofiler_output_bai = merlin_magic.tbprofiler_output_bai
    File? tbprofiler_output_vcf = merlin_magic.tbprofiler_output_vcf
    String? tbprofiler_version = merlin_magic.tbprofiler_version
    String? tbprofiler_main_lineage = merlin_magic.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = merlin_magic.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = merlin_magic.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = merlin_magic.tbprofiler_resistance_genes
    File? tbp_parser_lims_report_csv = merlin_magic.tbp_parser_lims_report_csv
    File? tbp_parser_looker_report_csv = merlin_magic.tbp_parser_looker_report_csv
    File? tbp_parser_laboratorian_report_csv = merlin_magic.tbp_parser_laboratorian_report_csv
    File? tbp_parser_coverage_report = merlin_magic.tbp_parser_coverage_report
    Float? tbp_parser_genome_percent_coverage = merlin_magic.tbp_parser_genome_percent_coverage
    Float? tbp_parser_average_genome_depth = merlin_magic.tbp_parser_average_genome_depth
    # Legionella pneumophila typing
    File? legsta_results = merlin_magic.legsta_results
    String? legsta_predicted_sbt = merlin_magic.legsta_predicted_sbt
    String? legsta_version = merlin_magic.legsta_version
    # Staphylococcus aureus
    File? spatyper_tsv = merlin_magic.spatyper_tsv
    String? spatyper_docker = merlin_magic.spatyper_docker
    String? spatyper_repeats = merlin_magic.spatyper_repeats
    String? spatyper_type = merlin_magic.spatyper_type
    String? spatyper_version = merlin_magic.spatyper_version
    File? staphopiasccmec_results_tsv = merlin_magic.staphopiasccmec_results_tsv
    File? staphopiasccmec_hamming_distance_tsv = merlin_magic.staphopiasccmec_hamming_distance_tsv
    String? staphopiasccmec_types_and_mecA_presence = merlin_magic.staphopiasccmec_types_and_mecA_presence
    String? staphopiasccmec_version = merlin_magic.staphopiasccmec_version
    String? staphopiasccmec_docker = merlin_magic.staphopiasccmec_docker
    File? agrvate_summary = merlin_magic.agrvate_summary
    File? agrvate_results = merlin_magic.agrvate_results
    String? agrvate_agr_group = merlin_magic.agrvate_agr_group
    String? agrvate_agr_match_score = merlin_magic.agrvate_agr_match_score
    String? agrvate_agr_canonical = merlin_magic.agrvate_agr_canonical
    String? agrvate_agr_multiple = merlin_magic.agrvate_agr_multiple
    String? agrvate_agr_num_frameshifts = merlin_magic.agrvate_agr_num_frameshifts
    String? agrvate_version = merlin_magic.agrvate_version
    String? agrvate_docker = merlin_magic.agrvate_docker
    # Streptococcus pneumoniae Typing
    String? pbptyper_predicted_1A_2B_2X = merlin_magic.pbptyper_predicted_1A_2B_2X
    File? pbptyper_pbptype_predicted_tsv = merlin_magic.pbptyper_pbptype_predicted_tsv
    String? pbptyper_version = merlin_magic.pbptyper_version
    String? pbptyper_docker = merlin_magic.pbptyper_docker
    String? poppunk_gps_cluster = merlin_magic.poppunk_gps_cluster
    File? poppunk_gps_external_cluster_csv = merlin_magic.poppunk_gps_external_cluster_csv
    String? poppunk_GPS_db_version = merlin_magic.poppunk_GPS_db_version
    String? poppunk_version = merlin_magic.poppunk_version
    String? poppunk_docker = merlin_magic.poppunk_docker
    # Streptococcus pyogenes Typing
    String? emmtyper_emm_type = merlin_magic.emmtyper_emm_type
    File? emmtyper_results_tsv = merlin_magic.emmtyper_results_tsv
    String? emmtyper_version = merlin_magic.emmtyper_version
    String? emmtyper_docker = merlin_magic.emmtyper_docker
    # Haemophilus influenzae Typing
    String? hicap_serotype = merlin_magic.hicap_serotype
    String? hicap_genes = merlin_magic.hicap_genes
    File? hicap_results_tsv = merlin_magic.hicap_results_tsv
    String? hicap_version = merlin_magic.hicap_version
    String? hicap_docker = merlin_magic.hicap_docker
    # Vibrio Typing
    File? srst2_vibrio_detailed_tsv = merlin_magic.srst2_vibrio_detailed_tsv
    String? srst2_vibrio_docker = merlin_magic.srst2_vibrio_docker
    String? srst2_vibrio_database = merlin_magic.srst2_vibrio_database
    String? srst2_vibrio_version = merlin_magic.srst2_vibrio_version
    String? srst2_vibrio_ctxA = merlin_magic.srst2_vibrio_ctxA
    String? srst2_vibrio_ompW = merlin_magic.srst2_vibrio_ompW
    String? srst2_vibrio_toxR = merlin_magic.srst2_vibrio_toxR
    String? srst2_vibrio_biotype = merlin_magic.srst2_vibrio_biotype
    String? srst2_vibrio_serogroup = merlin_magic.srst2_vibrio_serogroup
    File? abricate_vibrio_detailed_tsv = merlin_magic.abricate_vibrio_detailed_tsv
    String? abricate_vibrio_database = merlin_magic.abricate_vibrio_database
    String? abricate_vibrio_docker = merlin_magic.abricate_vibrio_docker
    String? abricate_vibrio_version = merlin_magic.abricate_vibrio_version
    String? abricate_vibrio_ctxA = merlin_magic.abricate_vibrio_ctxA
    String? abricate_vibrio_ompW = merlin_magic.abricate_vibrio_ompW
    String? abricate_vibrio_toxR = merlin_magic.abricate_vibrio_toxR
    String? abricate_vibrio_biotype = merlin_magic.abricate_vibrio_biotype
    String? abricate_vibrio_serogroup = merlin_magic.abricate_vibrio_serogroup
    # export taxon table output
    String? taxon_table_status = export_taxon_table.status
    # ARLN required outputs
    String? arln_r1_q30_raw = arln_stats.read1_raw_q30
    String? arln_r1_q30_clean = arln_stats.read1_clean_q30
    String? arln_assembly_ratio = arln_stats.assembly_ratio
    String? arln_taxon_assembly_ratio_stdev = arln_stats.taxon_assembly_ratio_stdev
    String? arln_taxon_gc_percent_stdev = arln_stats.taxon_gc_percent_stdev
    String? arln_taxon_gc_mean = arln_stats.taxon_gc_mean
    String? arln_assembly_zscore = arln_stats.assembly_zscore
    String? arln_stats_docker_version = arln_stats.docker_version
  }
}