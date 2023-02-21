version 1.0

import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc_workflow
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../../tasks/assembly/task_dragonflye.wdl" as dragonflye_task
import "../../tasks/quality_control/task_quast.wdl" as quast_task
import "../../tasks/quality_control/task_cg_pipeline.wdl" as cg_pipeline_task
import "../../tasks/quality_control/task_screen.wdl" as screen_task
import "../../tasks/quality_control/task_busco.wdl" as busco_task
import "../../tasks/taxon_id/task_gambit.wdl" as gambit_task
import "../../tasks/quality_control/task_mummer_ani.wdl" as ani_task
import "../../tasks/gene_typing/task_amrfinderplus.wdl" as amrfinderplus_task
import "../../tasks/gene_typing/task_resfinder.wdl" as resfinder_task
import "../../tasks/species_typing/task_ts_mlst.wdl" as ts_mlst_task
import "../../tasks/gene_typing/task_bakta.wdl" as bakta_task
import "../../tasks/gene_typing/task_prokka.wdl" as prokka_task
import "../../tasks/gene_typing/task_plasmidfinder.wdl" as plasmidfinder_task
import "../../tasks/quality_control/task_qc_check.wdl" as qc_check_task
import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/utilities/task_broad_terra_tools.wdl" as terra_tools_task

workflow theiaprok_ont {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of ONT bacterial NGS data"
  }
  input {
    String samplename
    String seq_method = "ONT"
    File read1
    Int? genome_size
    String? run_id
    String? collection_date
    String? originating_lab
    String? city
    String? county
    String? zip
    File? taxon_tables
    String terra_project = "NA"
    String terra_workspace = "NA"
    # by default do not call ANI task, but user has ability to enable this task if working with enteric pathogens or supply their own high-quality reference genome
    Boolean call_ani = false
    Int min_reads = 7472
    Int min_basepairs = 2241820
    Int min_genome_size = 100000
    Int max_genome_size = 18040666 
    Int min_coverage = 10
    Boolean call_resfinder = false
    Boolean skip_screen = false 
    String genome_annotation = "prokka"
    File? qc_check_table
    String? expected_taxon
  }
  call versioning_task.version_capture{
    input:
  }
  call screen_task.check_reads_ont as raw_check_reads {
    input:
      read1 = read1,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_size = min_genome_size,
      max_genome_size = max_genome_size,
      min_coverage = min_coverage,
      skip_screen = skip_screen
  }
  if (raw_check_reads.read_screen == "PASS") {
    call read_qc_workflow.read_QC_trim_ont as read_QC_trim {
      input:
        samplename = samplename,
        read1 = read1,
        genome_size = genome_size
    }
    call screen_task.check_reads_ont as clean_check_reads {
      input:
        read1 = read_QC_trim.read1_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_size = min_genome_size,
        max_genome_size = max_genome_size,
        min_coverage = min_coverage,
        skip_screen = skip_screen
    }
    if (clean_check_reads.read_screen == "PASS") {
       call dragonflye_task.dragonflye {
         input:
           read1 = read_QC_trim.read1_clean,
           genome_size = select_first([genome_size, read_QC_trim.est_genome_size]),
           samplename = samplename
       }
      call quast_task.quast {
        input:
          assembly = dragonflye.assembly_fasta,
          samplename = samplename
      }
      call cg_pipeline_task.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename,
          genome_length = select_first([genome_size, quast.genome_length])
      }
      call cg_pipeline_task.cg_pipeline as cg_pipeline_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename,
          genome_length = select_first([genome_size, quast.genome_length])
      }
      call gambit_task.gambit {
        input:
          assembly = dragonflye.assembly_fasta,
          samplename = samplename
      }
      call busco_task.busco {
        input:
          assembly = dragonflye.assembly_fasta,
          samplename = samplename
      }
      if (call_ani) {
        call ani_task.animummer as ani {
          input:
            assembly = dragonflye.assembly_fasta,
            samplename = samplename
        }
      }
      call amrfinderplus_task.amrfinderplus_nuc as amrfinderplus {
        input:
          assembly = dragonflye.assembly_fasta,
          samplename = samplename,
          organism = gambit.gambit_predicted_taxon
      }
      if (call_resfinder) {
        call resfinder_task.resfinder {
          input:
            assembly = dragonflye.assembly_fasta,
            samplename = samplename,
            organism = gambit.gambit_predicted_taxon
        }
      }
      call ts_mlst_task.ts_mlst {
        input: 
          assembly = dragonflye.assembly_fasta,
          samplename = samplename
      }
      if (genome_annotation == "prokka") {
        call prokka_task.prokka {
          input:
            assembly = dragonflye.assembly_fasta,
            samplename = samplename
        }
      }
      if (genome_annotation == "bakta") {
        call bakta_task.bakta {
          input:
            assembly = dragonflye.assembly_fasta,
            samplename = samplename
        }
      }
      call plasmidfinder_task.plasmidfinder {
        input:
          assembly = dragonflye.assembly_fasta,
          samplename = samplename
      }
      if(defined(qc_check_table)) {
        call qc_check_task.qc_check { # will request shelly's help in the future to make this applicable
          input:
            qc_check_table = qc_check_table,
            expected_taxon = expected_taxon,
            gambit_predicted_taxon = gambit.gambit_predicted_taxon,
            r1_mean_q_raw = cg_pipeline_raw.r1_mean_q,
            r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength,
            r1_mean_q_clean = cg_pipeline_clean.r1_mean_q,
            r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength,
            est_coverage_raw = cg_pipeline_raw.est_coverage,
            est_coverage_clean = cg_pipeline_clean.est_coverage,
            assembly_length = quast.genome_length,
            number_contigs = quast.number_contigs,
            n50_value = quast.n50_value,
            busco_results = busco.busco_results,
            ani_highest_percent = ani.ani_highest_percent,
            ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned,
            ani_top_species_match = ani.ani_top_species_match
        }
      }
      call merlin_magic_workflow.merlin_magic {
        input:
          merlin_tag = gambit.merlin_tag,
          assembly = dragonflye.assembly_fasta,
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          ont_data = true
      }
      if (defined(taxon_tables)) {
        call terra_tools_task.export_taxon_tables {
          input:
            terra_project = terra_project,
            terra_workspace = terra_workspace,
            sample_taxon = gambit.gambit_predicted_taxon,
            taxon_tables = taxon_tables,
            samplename = samplename,
            read1 = read1,
            read1_clean = read_QC_trim.read1_clean,
            run_id = run_id,
            collection_date = collection_date,
            originating_lab = originating_lab,
            city = city,
            county = county,
            zip = zip,
            theiaprok_ont_version = version_capture.phb_version,
            theiaprok_ont_analysis_date = version_capture.date,
            seq_platform = seq_method,
            num_reads_raw1 = read_QC_trim.number_raw_reads,
            fastq_scan_version = read_QC_trim.fastq_scan_version,
            num_reads_clean1 = read_QC_trim.number_clean_reads,
            r1_mean_q_raw = cg_pipeline_raw.r1_mean_q, 
            r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength,
            nanoq_version = read_QC_trim.nanoq_version,
            nanoplot_html = read_QC_trim.nanoplot_html,
            nanoplot_version = read_QC_trim.nanoplot_version,
            kmc_est_genome_size = read_QC_trim.est_genome_size,
            kmc_kmer_stats = read_QC_trim.kmc_kmer_stats,
            kmc_version = read_QC_trim.kmc_version,
            rasusa_version = read_QC_trim.rasusa_version,
            tiptoft_plasmid_replicon_fastq = read_QC_trim.tiptoft_plasmid_replicon_fastq,
            tiptoft_plasmid_replicon_genes = read_QC_trim.tiptoft_plasmid_replicon_genes,
            tiptoft_version = read_QC_trim.tiptoft_version,
            assembly_fasta = dragonflye.assembly_fasta,
            dragonflye_version = dragonflye.dragonflye_version,
            quast_report = quast.quast_report,
            quast_version = quast.version,
            assembly_length = quast.genome_length,
            number_contigs = quast.number_contigs,
            n50_value = quast.n50_value,
            quast_gc_percent = quast.gc_percent,
            cg_pipeline_report_raw = cg_pipeline_raw.cg_pipeline_report,
            cg_pipeline_docker = cg_pipeline_raw.cg_pipeline_docker,
            est_coverage_raw = cg_pipeline_raw.est_coverage,
            cg_pipeline_report_clean = cg_pipeline_clean.cg_pipeline_report,
            est_coverage_clean = cg_pipeline_clean.est_coverage,
            gambit_report = gambit.gambit_report_file,
            gambit_predicted_taxon = gambit.gambit_predicted_taxon,
            gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank,
            gambit_closest_genomes = gambit.gambit_closest_genomes_file,
            gambit_version = gambit.gambit_version,
            gambit_db_version = gambit.gambit_db_version,
            gambit_docker = gambit.gambit_docker,
            busco_version = busco.busco_version,
            busco_database = busco.busco_database,
            busco_results = busco.busco_results,
            busco_report = busco.busco_report,
            ani_highest_percent = ani.ani_highest_percent,
            ani_highest_percent_bases_aligned = ani.ani_highest_percent_bases_aligned,
            ani_output_tsv = ani.ani_output_tsv,
            ani_top_species_match = ani.ani_top_species_match,
            ani_mummer_version = ani.ani_mummer_version,
            amrfinderplus_all_report = amrfinderplus.amrfinderplus_all_report,
            amrfinderplus_amr_report = amrfinderplus.amrfinderplus_amr_report,
            amrfinderplus_stress_report = amrfinderplus.amrfinderplus_stress_report,
            amrfinderplus_virulence_report = amrfinderplus.amrfinderplus_virulence_report,
            amrfinderplus_amr_genes = amrfinderplus.amrfinderplus_amr_genes,
            amrfinderplus_stress_genes = amrfinderplus.amrfinderplus_stress_genes,
            amrfinderplus_virulence_genes = amrfinderplus.amrfinderplus_virulence_genes,
            amrfinderplus_amr_classes = amrfinderplus.amrfinderplus_amr_classes,
            amrfinderplus_amr_subclasses = amrfinderplus.amrfinderplus_amr_subclasses,
            amrfinderplus_version = amrfinderplus.amrfinderplus_version,
            amrfinderplus_db_version = amrfinderplus.amrfinderplus_db_version,
            resfinder_pheno_table = resfinder.resfinder_pheno_table,
            resfinder_pheno_table_species = resfinder.resfinder_pheno_table_species,
            resfinder_seqs = resfinder.resfinder_hit_in_genome_seq,
            resfinder_results = resfinder.resfinder_results_tab,
            resfinder_pointfinder_pheno_table = resfinder.pointfinder_pheno_table,
            resfinder_pointfinder_results = resfinder.pointfinder_results,
            resfinder_db_version = resfinder.resfinder_db_version,
            resfinder_docker = resfinder.resfinder_docker,
            ts_mlst_results = ts_mlst.ts_mlst_results,
            ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st,
            ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme,
            ts_mlst_version = ts_mlst.ts_mlst_version,
            ts_mlst_novel_alleles = ts_mlst.ts_mlst_novel_alleles,
            serotypefinder_report = merlin_magic.serotypefinder_report,
            serotypefinder_docker = merlin_magic.serotypefinder_docker,
            serotypefinder_serotype = merlin_magic.serotypefinder_serotype,
            ectyper_results = merlin_magic.ectyper_results,
            ectyper_version = merlin_magic.ectyper_version,
            ectyper_predicted_serotype = merlin_magic.ectyper_predicted_serotype,
            shigatyper_predicted_serotype = merlin_magic.shigatyper_predicted_serotype,
            shigatyper_ipaB_presence_absence = merlin_magic.shigatyper_ipaB_presence_absence,
            shigatyper_notes = merlin_magic.shigatyper_notes,
            shigatyper_hits_tsv = merlin_magic.shigatyper_hits_tsv,
            shigatyper_summary_tsv = merlin_magic.shigatyper_summary_tsv,
            shigatyper_version = merlin_magic.shigatyper_version,
            shigatyper_docker = merlin_magic.shigatyper_docker,
            shigeifinder_report = merlin_magic.shigeifinder_report,
            shigeifinder_docker = merlin_magic.shigeifinder_docker,
            shigeifinder_version = merlin_magic.shigeifinder_version,
            shigeifinder_ipaH_presence_absence = merlin_magic.shigeifinder_ipaH_presence_absence,
            shigeifinder_num_virulence_plasmid_genes = merlin_magic.shigeifinder_num_virulence_plasmid_genes,
            shigeifinder_cluster = merlin_magic.shigeifinder_cluster,
            shigeifinder_serotype = merlin_magic.shigeifinder_serotype,
            shigeifinder_O_antigen = merlin_magic.shigeifinder_O_antigen,
            shigeifinder_H_antigen = merlin_magic.shigeifinder_H_antigen,
            shigeifinder_notes = merlin_magic.shigeifinder_notes,
            sonneityping_mykrobe_report_csv = merlin_magic.sonneityping_mykrobe_report_csv,
            sonneityping_mykrobe_report_json = merlin_magic.sonneityping_mykrobe_report_json,
            sonneityping_final_report_tsv = merlin_magic.sonneityping_final_report_tsv,
            sonneityping_mykrobe_version = merlin_magic.sonneityping_mykrobe_version,
            sonneityping_mykrobe_docker = merlin_magic.sonneityping_mykrobe_docker,
            sonneityping_species = merlin_magic.sonneityping_species,
            sonneityping_final_genotype = merlin_magic.sonneityping_final_genotype,
            sonneityping_genotype_confidence = merlin_magic.sonneityping_genotype_confidence,
            sonneityping_genotype_name = merlin_magic.sonneityping_genotype_name,
            lissero_results = merlin_magic.lissero_results,
            lissero_version = merlin_magic.lissero_version,
            lissero_serotype = merlin_magic.lissero_serotype,
            sistr_results = merlin_magic.sistr_results,
            sistr_allele_json = merlin_magic.sistr_allele_json,
            sister_allele_fasta = merlin_magic.sistr_allele_fasta,
            sistr_cgmlst = merlin_magic.sistr_cgmlst,
            sistr_version = merlin_magic.sistr_version,
            sistr_predicted_serotype = merlin_magic.sistr_predicted_serotype,
            seqsero2_report = merlin_magic.seqsero2_report,
            seqsero2_version = merlin_magic.seqsero2_version,
            seqsero2_predicted_antigenic_profile = merlin_magic.seqsero2_predicted_antigenic_profile,
            seqsero2_predicted_serotype = merlin_magic.seqsero2_predicted_serotype,
            seqsero2_predicted_contamination = merlin_magic.seqsero2_predicted_contamination,
            genotyphi_report_tsv = merlin_magic.genotyphi_report_tsv,
            genotyphi_mykrobe_json = merlin_magic.genotyphi_mykrobe_json,
            genotyphi_version = merlin_magic.genotyphi_version,
            genotyphi_species = merlin_magic.genotyphi_species,
            genotyphi_st_probes_percent_coverage = merlin_magic.genotyphi_st_probes_percent_coverage,
            genotyphi_final_genotype = merlin_magic.genotyphi_final_genotype,
            genotyphi_genotype_confidence = merlin_magic.genotyphi_genotype_confidence,
            kleborate_output_file = merlin_magic.kleborate_output_file,
            kleborate_version = merlin_magic.kleborate_version,
            kleborate_docker = merlin_magic.kleborate_docker,
            kleborate_key_resistance_genes = merlin_magic.kleborate_key_resistance_genes,
            kleborate_genomic_resistance_mutations = merlin_magic.kleborate_genomic_resistance_mutations,
            kleborate_mlst_sequence_type = merlin_magic.kleborate_mlst_sequence_type,
            kleborate_klocus = merlin_magic.kleborate_klocus,
            kleborate_ktype = merlin_magic.kleborate_ktype,
            kleborate_olocus = merlin_magic.kleborate_olocus,
            kleborate_otype = merlin_magic.kleborate_otype,
            kleborate_klocus_confidence = merlin_magic.kleborate_klocus_confidence,
            kleborate_olocus_confidence = merlin_magic.kleborate_olocus_confidence,
            kaptive_output_file_k = merlin_magic.kaptive_output_file_k,
            kaptive_output_file_oc = merlin_magic.kaptive_output_file_oc,
            kaptive_version = merlin_magic.kaptive_version,
            kaptive_k_locus = merlin_magic.kaptive_k_match,
            kaptive_k_type = merlin_magic.kaptive_k_type,
            kaptive_kl_confidence = merlin_magic.kaptive_k_confidence,
            kaptive_oc_locus = merlin_magic.kaptive_oc_match,
            kaptive_ocl_confidence = merlin_magic.kaptive_oc_confidence,
            abricate_abaum_plasmid_tsv = merlin_magic.abricate_results,
            abricate_abaum_plasmid_type_genes = merlin_magic.abricate_genes,
            abricate_database = merlin_magic.abricate_database,
            abricate_version = merlin_magic.abricate_version,
            abricate_docker = merlin_magic.abricate_docker,
            tbprofiler_output_file = merlin_magic.tbprofiler_output_file,
            tbprofiler_output_bam = merlin_magic.tbprofiler_output_bam,
            tbprofiler_output_bai = merlin_magic.tbprofiler_output_bai,
            tbprofiler_version = merlin_magic.tbprofiler_version,
            tbprofiler_main_lineage = merlin_magic.tbprofiler_main_lineage,
            tbprofiler_sub_lineage = merlin_magic.tbprofiler_sub_lineage,
            tbprofiler_dr_type = merlin_magic.tbprofiler_dr_type,
            tbprofiler_resistance_genes = merlin_magic.tbprofiler_resistance_genes,
            legsta_results = merlin_magic.legsta_results,
            legsta_predicted_sbt = merlin_magic.legsta_predicted_sbt,
            legsta_version = merlin_magic.legsta_version,
            prokka_gff = prokka.prokka_gff,
            prokka_gbk = prokka.prokka_gbk,
            prokka_sqn = prokka.prokka_sqn,
            bakta_gbff = bakta.bakta_gbff,
            bakta_gff3 = bakta.bakta_gff3,
            bakta_tsv = bakta.bakta_tsv,
            bakta_summary = bakta.bakta_txt,
            bakta_version = bakta.bakta_version,
            plasmidfinder_plasmids = plasmidfinder.plasmidfinder_plasmids,
            plasmidfinder_results = plasmidfinder.plasmidfinder_results,
            plasmidfinder_seqs = plasmidfinder.plasmidfinder_seqs,
            plasmidfinder_docker = plasmidfinder.plasmidfinder_docker,
            plasmidfinder_db_version = plasmidfinder.plasmidfinder_db_version,
            pbptyper_predicted_1A_2B_2X = merlin_magic.pbptyper_predicted_1A_2B_2X,
            pbptyper_pbptype_predicted_tsv = merlin_magic.pbptyper_pbptype_predicted_tsv,
            pbptyper_version = merlin_magic.pbptyper_version,
            pbptyper_docker = merlin_magic.pbptyper_docker,
            poppunk_gps_cluster = merlin_magic.poppunk_gps_cluster,
            poppunk_gps_external_cluster_csv = merlin_magic.poppunk_gps_external_cluster_csv,
            poppunk_GPS_db_version = merlin_magic.poppunk_GPS_db_version,
            poppunk_version = merlin_magic.poppunk_version,
            poppunk_docker = merlin_magic.poppunk_docker,
            pasty_serogroup = merlin_magic.pasty_serogroup,
            pasty_serogroup_coverage = merlin_magic.pasty_serogroup_coverage,
            pasty_serogroup_fragments = merlin_magic.pasty_serogroup_fragments,
            pasty_summary_tsv = merlin_magic.pasty_summary_tsv,
            pasty_blast_hits = merlin_magic.pasty_blast_hits,
            pasty_all_serogroups = merlin_magic.pasty_all_serogroups,
            pasty_version = merlin_magic.pasty_version,
            pasty_docker = merlin_magic.pasty_docker,
            pasty_comment = merlin_magic.pasty_comment,
            qc_check = qc_check.qc_check,
            qc_standard = qc_check.qc_standard
        }
      }
    }
  }
  output {
    # Version Captures
    String theiaprok_ont_version = version_capture.phb_version
    String theiaprok_ont_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Read QC - fastq_scan and nanoq outputs
    File? read1_clean = read_QC_trim.read1_clean
    Int? num_reads_raw1 = read_QC_trim.number_raw_reads
    Int? num_reads_clean1 = read_QC_trim.number_clean_reads
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    String? nanoq_version = read_QC_trim.nanoq_version
    # Read QC - nanoplot outputs
    File? nanoplot_html = read_QC_trim.nanoplot_html
    String? nanoplot_version = read_QC_trim.nanoplot_version
    # Read QC - kmc outputs
    String? kmc_est_genome_size = read_QC_trim.est_genome_size
    File? kmc_kmer_stats = read_QC_trim.kmc_kmer_stats
    String? kmc_version = read_QC_trim.kmc_version
    # Read QC - rasusa outputs
    String? rasusa_version = read_QC_trim.rasusa_version
    # Read QC - tiptoft outputs
    File? tiptoft_plasmid_replicon_fastq = read_QC_trim.tiptoft_plasmid_replicon_fastq
    String? tiptoft_plasmid_replicon_genes = read_QC_trim.tiptoft_plasmid_replicon_genes
    String? tiptoft_version = read_QC_trim.tiptoft_version
    # Read QC - cg pipeline outputs
    Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
    Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
    # Assembly - dragonflye outputs
    File? assembly_fasta = dragonflye.assembly_fasta
    String? dragonflye_version = dragonflye.dragonflye_version
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
    # NCBI-AMRFinderPlus Outputs
    File? amrfinderplus_all_report = amrfinderplus.amrfinderplus_all_report
    File? amrfinderplus_amr_report = amrfinderplus.amrfinderplus_amr_report
    File? amrfinderplus_stress_report = amrfinderplus.amrfinderplus_stress_report
    File? amrfinderplus_virulence_report = amrfinderplus.amrfinderplus_virulence_report
    String? amrfinderplus_amr_genes = amrfinderplus.amrfinderplus_amr_genes
    String? amrfinderplus_stress_genes = amrfinderplus.amrfinderplus_stress_genes
    String? amrfinderplus_virulence_genes = amrfinderplus.amrfinderplus_virulence_genes
    String? amrfinderplus_amr_classes = amrfinderplus.amrfinderplus_amr_classes
    String? amrfinderplus_amr_subclasses = amrfinderplus.amrfinderplus_amr_subclasses
    String? amrfinderplus_version = amrfinderplus.amrfinderplus_version
    String? amrfinderplus_db_version = amrfinderplus.amrfinderplus_db_version
    # Resfinder Outputs
    File? resfinder_pheno_table = resfinder.resfinder_pheno_table
    File? resfinder_pheno_table_species = resfinder.resfinder_pheno_table_species
    File? resfinder_seqs = resfinder.resfinder_hit_in_genome_seq
    File? resfinder_results = resfinder.resfinder_results_tab
    File? resfinder_pointfinder_pheno_table = resfinder.pointfinder_pheno_table
    File? resfinder_pointfinder_results = resfinder.pointfinder_results
    String? resfinder_db_version = resfinder.resfinder_db_version
    String? resfinder_docker = resfinder.resfinder_docker
    # MLST Typing
    File? ts_mlst_results = ts_mlst.ts_mlst_results
    String? ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st
    String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
    String? ts_mlst_version = ts_mlst.ts_mlst_version
    File? ts_mlst_novel_alleles = ts_mlst.ts_mlst_novel_alleles
    # Prokka Results
    File? prokka_gff = prokka.prokka_gff
    File? prokka_gbk = prokka.prokka_gbk
    File? prokka_sqn = prokka.prokka_sqn
    # Bakta Results
    File? bakta_gbff = bakta.bakta_gbff
    File? bakta_gff3 = bakta.bakta_gff3
    File? bakta_tsv = bakta.bakta_tsv
    File? bakta_summary = bakta.bakta_txt
    String? bakta_version = bakta.bakta_version
    # Plasmidfinder Results
    String? plasmidfinder_plasmids = plasmidfinder.plasmidfinder_plasmids
    File? plasmidfinder_results = plasmidfinder.plasmidfinder_results
    File? plasmidfinder_seqs = plasmidfinder.plasmidfinder_seqs
    String? plasmidfinder_docker = plasmidfinder.plasmidfinder_docker
    String? plasmidfinder_db_version = plasmidfinder.plasmidfinder_db_version
    # Ecoli Typing
    File? serotypefinder_report = merlin_magic.serotypefinder_report
    String? serotypefinder_docker = merlin_magic.serotypefinder_docker
    String? serotypefinder_serotype = merlin_magic.serotypefinder_serotype
    File? ectyper_results = merlin_magic.ectyper_results
    String? ectyper_version = merlin_magic.ectyper_version
    String? ectyper_predicted_serotype = merlin_magic.ectyper_predicted_serotype
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
    File? sister_allele_fasta = merlin_magic.sistr_allele_fasta
    File? sistr_cgmlst = merlin_magic.sistr_cgmlst
    String? sistr_version = merlin_magic.sistr_version
    String? sistr_predicted_serotype = merlin_magic.sistr_predicted_serotype
    File? seqsero2_report = merlin_magic.seqsero2_report
    String? seqsero2_version = merlin_magic.seqsero2_version
    String? seqsero2_predicted_antigenic_profile = merlin_magic.seqsero2_predicted_antigenic_profile
    String? seqsero2_predicted_serotype = merlin_magic.seqsero2_predicted_serotype
    String? seqsero2_predicted_contamination = merlin_magic.seqsero2_predicted_contamination
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
    # Acinetobacter Typing
    File? kaptive_output_file_k = merlin_magic.kaptive_output_file_k
    File? kaptive_output_file_oc = merlin_magic.kaptive_output_file_oc
    String? kaptive_version = merlin_magic.kaptive_version
    String? kaptive_k_locus = merlin_magic.kaptive_k_match
    String? kaptive_k_type = merlin_magic.kaptive_k_type
    String? kaptive_kl_confidence = merlin_magic.kaptive_k_confidence
    String? kaptive_oc_locus = merlin_magic.kaptive_oc_match
    String? kaptive_ocl_confidence = merlin_magic.kaptive_oc_confidence
    File? abricate_abaum_plasmid_tsv = merlin_magic.abricate_results
    String? abricate_abaum_plasmid_type_genes = merlin_magic.abricate_genes
    String? abricate_database = merlin_magic.abricate_database
    String? abricate_version = merlin_magic.abricate_version
    String? abricate_docker = merlin_magic.abricate_docker
    # Mycobacterium Typing
    File? tbprofiler_output_file = merlin_magic.tbprofiler_output_file
    File? tbprofiler_output_bam = merlin_magic.tbprofiler_output_bam
    File? tbprofiler_output_bai = merlin_magic.tbprofiler_output_bai
    String? tbprofiler_version = merlin_magic.tbprofiler_version
    String? tbprofiler_main_lineage = merlin_magic.tbprofiler_main_lineage
    String? tbprofiler_sub_lineage = merlin_magic.tbprofiler_sub_lineage
    String? tbprofiler_dr_type = merlin_magic.tbprofiler_dr_type
    String? tbprofiler_resistance_genes = merlin_magic.tbprofiler_resistance_genes
    # Legionella pneumophila typing
    File? legsta_results = merlin_magic.legsta_results
    String? legsta_predicted_sbt = merlin_magic.legsta_predicted_sbt
    String? legsta_version = merlin_magic.legsta_version
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
    # export taxon table output
    String? taxon_table_status = export_taxon_tables.status
  }
}