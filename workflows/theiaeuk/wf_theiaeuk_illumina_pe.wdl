version 1.0

import "../../workflows/utilities/wf_digger_denovo.wdl" as digger_denovo
import "../../tasks/quality_control/advanced_metrics/task_busco.wdl" as busco_task
import "../../tasks/quality_control/basic_statistics/task_cg_pipeline.wdl" as cg_pipeline_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_gambit.wdl" as gambit_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../utilities/wf_read_QC_trim_pe.wdl" as read_qc

workflow theiaeuk_illumina_pe {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of paired-end eukaryotic NGS data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1
    File read2
    Boolean call_rasusa = true
    Int min_reads = 30000
    # Edited default values
    Int min_basepairs = 45000000
    Int min_genome_length = 9000000
    Int max_genome_length = 178000000
    Int min_coverage = 10
    Int min_proportion = 40
    Int trim_min_length = 75
    Int trim_quality_min_score = 20
    Int trim_window_size = 10
    Int min_contig_length = 1000
    Int busco_memory = 24
    String busco_docker_image = "us-docker.pkg.dev/general-theiagen/ezlabgva/busco:v5.3.2_cv1"
    Boolean skip_screen = false 
    File? qc_check_table
    String? expected_taxon
    Int? genome_length 
    Float subsample_coverage = 150 # default coverage for RASUSA is set to 150X
    Int cpu = 8
    Int memory = 16
    # default gambit outputs
    File gambit_db_genomes = "gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-metadata-1.0.0-20241213.gdb"
    File gambit_db_signatures = "gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-signatures-1.0.0-20241213.gs"
  }
  call versioning.version_capture {
    input:
  } 
  if (! skip_screen) {
    call screen.check_reads as raw_check_reads {
      input:
        read1 = read1,
        read2 = read2,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        workflow_series = "theiaeuk",
        expected_genome_length = genome_length
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    if (call_rasusa) {
      call rasusa.rasusa as rasusa_task {
        input:
          read1 = read1,
          read2 = read2,
          samplename = samplename,
          genome_length = select_first([genome_length, raw_check_reads.est_genome_length, 0]),
          coverage = subsample_coverage
      }
    }  
    call read_qc.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1 = select_first([rasusa_task.read1_subsampled, read1]),
        read2 = select_first([rasusa_task.read2_subsampled, read2]),
        trim_min_length = trim_min_length,
        trim_quality_min_score = trim_quality_min_score,
        trim_window_size = trim_window_size,
        workflow_series = "theiaeuk"
    }
    if (! skip_screen) {
      call screen.check_reads as clean_check_reads {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          min_reads = min_reads,
          min_basepairs = min_basepairs,
          min_genome_length = min_genome_length,
          max_genome_length = max_genome_length,
          min_coverage = min_coverage,
          min_proportion = min_proportion,
          workflow_series = "theiaeuk",
          expected_genome_length = genome_length
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      call digger_denovo.digger_denovo {
        input:
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          min_contig_length = min_contig_length,
          filter_contigs_min_length = min_contig_length
      }
      call quast_task.quast {
        input:
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename,
          cpu = cpu,
          memory = memory
      }
      call cg_pipeline_task.cg_pipeline as cg_pipeline_raw {
        input:
          read1 = read1,
          read2 = read2,
          samplename = samplename,
          genome_length = select_first([quast.genome_length, clean_check_reads.est_genome_length]),
          cpu = cpu,
          memory = memory
      }
      call cg_pipeline_task.cg_pipeline as cg_pipeline_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          samplename = samplename,
          genome_length = select_first([quast.genome_length, clean_check_reads.est_genome_length]),
          cpu = cpu,
          memory = memory
      }
      call gambit_task.gambit {
        input:
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename,
          gambit_db_genomes = gambit_db_genomes,
          gambit_db_signatures = gambit_db_signatures,
          cpu = cpu,
          memory = memory
      }
      call busco_task.busco {
        input:
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename,
          eukaryote = true,
          memory = busco_memory,
          docker = busco_docker_image
      }
      if (defined(qc_check_table)) {
        call qc_check.qc_check_phb as qc_check_task {
          input:
            qc_check_table = qc_check_table,
            expected_taxon = expected_taxon,
            gambit_predicted_taxon = gambit.gambit_predicted_taxon,
            num_reads_raw1 = read_QC_trim.fastq_scan_raw1,
            num_reads_raw2 = read_QC_trim.fastq_scan_raw2,
            num_reads_clean1 = read_QC_trim.fastq_scan_clean1,
            num_reads_clean2 = read_QC_trim.fastq_scan_clean2,
            r1_mean_q_raw = cg_pipeline_raw.r1_mean_q,
            r2_mean_q_raw = cg_pipeline_raw.r2_mean_q,
            combined_mean_q_raw = cg_pipeline_raw.combined_mean_q,
            r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength,
            r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength,  
            combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength,
            r1_mean_q_clean = cg_pipeline_clean.r1_mean_q,
            r2_mean_q_clean = cg_pipeline_clean.r2_mean_q,
            combined_mean_q_clean = cg_pipeline_clean.combined_mean_q,
            r1_mean_readlength_clean = cg_pipeline_clean.r1_mean_readlength,
            r2_mean_readlength_clean = cg_pipeline_clean.r2_mean_readlength,  
            combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength,    
            est_coverage_raw = cg_pipeline_raw.est_coverage,
            est_coverage_clean = cg_pipeline_clean.est_coverage,
            assembly_length = quast.genome_length,
            number_contigs = quast.number_contigs,
            n50_value = quast.n50_value,
            quast_gc_percent = quast.gc_percent,
            busco_results = busco.busco_results
        }
      }
      call merlin_magic_workflow.merlin_magic {
        input:
          merlin_tag = gambit.merlin_tag,
          assembly = digger_denovo.assembly_fasta,
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean,
          theiaeuk = true
      }
    }
  }
  output {
    # Version Captures
    String theiaeuk_illumina_pe_version = version_capture.phb_version
    String theiaeuk_illumina_pe_analysis_date = version_capture.date
    # RASUSA
    String? rasusa_version = rasusa_task.rasusa_version
    File? read1_subsampled = rasusa_task.read1_subsampled
    File? read2_subsampled = rasusa_task.read2_subsampled
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC - fastq_scan outputs
    Int? fastq_scan_num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? fastq_scan_num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? fastq_scan_num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? fastq_scan_num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? fastq_scan_num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? fastq_scan_num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    File? fastq_scan_raw1_json = read_QC_trim.fastq_scan_raw1_json
    File? fastq_scan_raw2_json = read_QC_trim.fastq_scan_raw2_json
    File? fastq_scan_clean1_json = read_QC_trim.fastq_scan_clean1_json
    File? fastq_scan_clean2_json = read_QC_trim.fastq_scan_clean2_json
    # Read QC - trimmomatic outputs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? trimmomatic_docker = read_QC_trim.trimmomatic_docker
    # Read QC - fastqc outputs
    Int? fastqc_num_reads_raw1 = read_QC_trim.fastqc_raw1
    Int? fastqc_num_reads_raw2 = read_QC_trim.fastqc_raw2
    String? fastqc_num_reads_raw_pairs = read_QC_trim.fastqc_raw_pairs
    Int? fastqc_num_reads_clean1 = read_QC_trim.fastqc_clean1
    Int? fastqc_num_reads_clean2 = read_QC_trim.fastqc_clean2
    String? fastqc_num_reads_clean_pairs = read_QC_trim.fastqc_clean_pairs
    File? fastqc_raw1_html = read_QC_trim.fastqc_raw1_html
    File? fastqc_raw2_html = read_QC_trim.fastqc_raw2_html
    File? fastqc_clean1_html = read_QC_trim.fastqc_clean1_html
    File? fastqc_clean2_html = read_QC_trim.fastqc_clean2_html
    String? fastqc_version = read_QC_trim.fastqc_version
    String? fastqc_docker = read_QC_trim.fastqc_docker
    # Read QC - fastp outputs
    String? fastp_version = read_QC_trim.fastp_version
    File? fastp_html_report = read_QC_trim.fastp_html_report
    # Read QC - bbduk outputs
    String? bbduk_docker = read_QC_trim.bbduk_docker
    File? read1_clean = read_QC_trim.read1_clean
    File? read2_clean = read_QC_trim.read2_clean
    # Read QC - cg pipeline outputs
    Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
    Float? r2_mean_q_raw = cg_pipeline_raw.r2_mean_q
    Float? combined_mean_q_raw = cg_pipeline_raw.combined_mean_q
    Float? combined_mean_q_clean = cg_pipeline_clean.combined_mean_q
    Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
    Float? r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength
    Float? combined_mean_readlength_raw = cg_pipeline_raw.combined_mean_readlength
    Float? combined_mean_readlength_clean = cg_pipeline_clean.combined_mean_readlength
    # Read QC - kraken outputs
    String? kraken2_version = read_QC_trim.kraken_version
    String? kraken2_report = read_QC_trim.kraken_report
    String? kraken2_database = read_QC_trim.kraken_database
    String? kraken_docker = read_QC_trim.kraken_docker
    # Assembly - shovill outputs and Assembly QC
    File? assembly_fasta = shovill_pe.assembly_fasta
    File? contigs_gfa = shovill_pe.contigs_gfa
    File? contigs_fastg = shovill_pe.contigs_fastg
    File? contigs_lastgraph = shovill_pe.contigs_lastgraph
    String? shovill_pe_version = shovill_pe.shovill_version
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
    # Taxon ID
    File? gambit_report = gambit.gambit_report_file
    File? gambit_closest_genomes = gambit.gambit_closest_genomes_file
    String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String? gambit_version = gambit.gambit_version
    String? gambit_db_version = gambit.gambit_db_version
    String? gambit_docker = gambit.gambit_docker
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # AMR_Search
    File? amr_search_results = merlin_magic.amr_search_results
    File? amr_search_csv = merlin_magic.amr_results_csv
    File? amr_search_results_pdf = merlin_magic.amr_results_pdf
    String? amr_search_docker = merlin_magic.amr_search_docker
    String? amr_search_version = merlin_magic.amr_search_version    
    # Cladetyper Outputs
    String? cladetyper_clade = merlin_magic.clade_type
    String? cladetyper_gambit_version = merlin_magic.cladetyper_version
    String? cladetyper_docker_image = merlin_magic.cladetyper_docker_image
    String? cladetyper_annotated_reference = merlin_magic.cladetype_annotated_ref
    # Snippy Outputs
    String? theiaeuk_snippy_variants_version = merlin_magic.snippy_variants_version
    String? theiaeuk_snippy_variants_query = merlin_magic.snippy_variants_query
    String? theiaeuk_snippy_variants_query_check = merlin_magic.snippy_variants_query_check
    String? theiaeuk_snippy_variants_hits = merlin_magic.snippy_variants_hits
    String? theiaeuk_snippy_variants_reference_genome = merlin_magic.snippy_variants_reference_genome
    String? theiaeuk_snippy_variants_gene_query_results = merlin_magic.snippy_variants_gene_query_results
    # Array[File]? snippy_outputs = merlin_magic.snippy_outputs
    String? theiaeuk_snippy_variants_results = merlin_magic.snippy_variants_results
    String? theiaeuk_snippy_variants_bam = merlin_magic.snippy_variants_bam
    String? theiaeuk_snippy_variants_bai = merlin_magic.snippy_variants_bai
    String? theiaeuk_snippy_variants_outdir_tarball = merlin_magic.snippy_variants_outdir_tarball
    String? theiaeuk_snippy_variants_summary = merlin_magic.snippy_variants_summary
    String? theiaeuk_snippy_variants_num_reads_aligned = merlin_magic.snippy_variants_num_reads_aligned
    String? theiaeuk_snippy_variants_coverage_tsv = merlin_magic.snippy_variants_coverage_tsv
    String? theiaeuk_snippy_variants_num_variants = merlin_magic.snippy_variants_num_variants
    String? theiaeuk_snippy_variants_percent_ref_coverage = merlin_magic.snippy_variants_percent_ref_coverage
  }
}
