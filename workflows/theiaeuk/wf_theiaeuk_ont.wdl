version 1.0

import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc
import "../utilities/wf_flye_denovo.wdl" as flye_workflow
import "../../tasks/taxon_id/task_gambit.wdl" as gambit
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/advanced_metrics/task_busco.wdl" as busco_task
import "../../tasks/quality_control/comparisons/task_screen.wdl" as screen_task

workflow theiaeuk_ont {
  input {
    File read1
    String samplename
    Int genome_length = 50000000 
    String workflow_series = "theiaeuk"
    Int busco_memory = 24
    String busco_docker_image = "us-docker.pkg.dev/general-theiagen/ezlabgva/busco:v5.3.2_cv1"
    File gambit_db_genomes = "gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-metadata-1.0.0-20241213.gdb"
    File gambit_db_signatures = "gs://gambit-databases-rp/fungal-version/1.0.0/gambit-fungal-signatures-1.0.0-20241213.gs"
    # read screen parameters
    Int min_reads = 5000
    Boolean skip_screen = false
    Boolean skip_mash = true
    Int min_basepairs = 45000000
    Int min_genome_length = 9000000
    Int max_genome_length = 178000000
    Int min_coverage = 5
  }
  call versioning.version_capture {
    input:
  }
  if (! skip_screen) {
    call screen_task.check_reads_se as raw_check_reads {
      input:
        read1 = read1,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_length = min_genome_length,
        max_genome_length = max_genome_length,
        min_coverage = min_coverage,
        skip_mash = skip_mash,
        expected_genome_length = genome_length
    }
  }
  if (select_first([raw_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
    call read_qc.read_QC_trim_ont as read_QC_trim {
      input:
        read1 = read1,
        samplename = samplename,
        genome_length = genome_length,
        workflow_series = workflow_series
    }
    if (! skip_screen) {
      call screen_task.check_reads_se as clean_check_reads {
        input:
          read1 = read_QC_trim.read1_clean,
          min_reads = min_reads,
          min_basepairs = min_basepairs,
          min_genome_length = min_genome_length,
          max_genome_length = max_genome_length,
          min_coverage = min_coverage,
          skip_mash = skip_mash,
          expected_genome_length = genome_length
      }
    }
    if (select_first([clean_check_reads.read_screen, ""]) == "PASS" || skip_screen) {
      call flye_workflow.flye_denovo {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename
      }
      #call quast on the assembly
      call quast_task.quast {
        input:
          assembly = flye_denovo.assembly_fasta,
          samplename = samplename
      }
      # nanoplot for raw reads
      call nanoplot_task.nanoplot as nanoplot_raw {
        input:
          read1 = read1,
          samplename = samplename,
          est_genome_length = select_first([genome_length, quast.genome_length])
      }
      # nanoplot for cleaned reads
      call nanoplot_task.nanoplot as nanoplot_clean {
        input:
          read1 = read_QC_trim.read1_clean,
          samplename = samplename,
          est_genome_length = select_first([genome_length, quast.genome_length])
      }
      # busco on the assembly
      call busco_task.busco {
        input:
          assembly = flye_denovo.assembly_fasta,
          samplename = samplename,
          eukaryote = true,
          memory = busco_memory,
          docker = busco_docker_image
      }
      # call gambit to predict taxon
      call gambit.gambit {
        input:
          assembly = flye_denovo.assembly_fasta,
          samplename = samplename,
          gambit_db_genomes = gambit_db_genomes,
          gambit_db_signatures = gambit_db_signatures
      }
      # call merlin magic for cladetyper and AMR search, snippy variants
       call merlin_magic_workflow.merlin_magic {
        input:
          samplename = samplename,
          merlin_tag = gambit.merlin_tag,
          assembly = flye_denovo.assembly_fasta,
          assembly_only = true,
          ont_data = true,
          theiaeuk = true
      }
    }
  }
  output {
    # Version Capture
    String theiaeuk_ont_version = version_capture.phb_version
    String theiaeuk_ont_analysis_date = version_capture.date
    # Sample Screening
    String? read_screen_raw = raw_check_reads.read_screen
    File? read_screen_raw_tsv = raw_check_reads.read_screen_tsv
    String? read_screen_clean = clean_check_reads.read_screen
    File? read_screen_clean_tsv = clean_check_reads.read_screen_tsv
    # Read QC outputs
    File? read1_clean = read_QC_trim.read1_clean
    String? nanoq_version = read_QC_trim.nanoq_version
    Int? est_genome_length = read_QC_trim.est_genome_length
    # Assembly - flye_denovo outputs
    File? assembly_fasta = flye_denovo.assembly_fasta
    File? contigs_gfa = flye_denovo.contigs_gfa
    File? bandage_plot = flye_denovo.bandage_plot
    File? filtered_contigs_metrics = flye_denovo.filtered_contigs_metrics
    String? flye_assembly_info = flye_denovo.flye_assembly_info
    String? medaka_model = flye_denovo.medaka_model_used
    String? porechop_version = flye_denovo.porechop_version
    String? flye_version = flye_denovo.flye_version
    String? bandage_version = flye_denovo.bandage_version
    String? medaka_version = flye_denovo.medaka_version
    String? racon_version = flye_denovo.racon_version
    String? bwa_version = flye_denovo.bwa_version
    String? polypolish_version = flye_denovo.polypolish_version
    String? dnaapler_version = flye_denovo.dnaapler_version
    # Read QC - nanoplot raw outputs
    File? nanoplot_html_raw = nanoplot_raw.nanoplot_html
    File? nanoplot_tsv_raw = nanoplot_raw.nanoplot_tsv
    Int? nanoplot_num_reads_raw1 = nanoplot_raw.num_reads
    Float? nanoplot_r1_median_readlength_raw = nanoplot_raw.median_readlength
    Float? nanoplot_r1_mean_readlength_raw = nanoplot_raw.mean_readlength
    Float? nanoplot_r1_stdev_readlength_raw = nanoplot_raw.stdev_readlength
    Float? nanoplot_r1_n50_raw = nanoplot_raw.n50
    Float? nanoplot_r1_mean_q_raw = nanoplot_raw.mean_q
    Float? nanoplot_r1_median_q_raw = nanoplot_raw.median_q
    Float? nanoplot_r1_est_coverage_raw = nanoplot_raw.est_coverage
    # Read QC - nanoplot clean outputs
    File? nanoplot_html_clean = nanoplot_clean.nanoplot_html
    File? nanoplot_tsv_clean = nanoplot_clean.nanoplot_tsv
    Int? nanoplot_num_reads_clean1 = nanoplot_clean.num_reads
    Float? nanoplot_r1_median_readlength_clean = nanoplot_clean.median_readlength
    Float? nanoplot_r1_mean_readlength_clean = nanoplot_clean.mean_readlength
    Float? nanoplot_r1_stdev_readlength_clean = nanoplot_clean.stdev_readlength
    Float? nanoplot_r1_n50_clean = nanoplot_clean.n50
    Float? nanoplot_r1_mean_q_clean = nanoplot_clean.mean_q
    Float? nanoplot_r1_median_q_clean = nanoplot_clean.median_q
    Float? nanoplot_r1_est_coverage_clean = nanoplot_clean.est_coverage
    # Read QC - nanoplot general outputs
    String? nanoplot_version = nanoplot_raw.nanoplot_version
    String? nanoplot_docker = nanoplot_raw.nanoplot_docker
    # Assembly QC - quast outputs
    File? quast_report = quast.quast_report
    String? quast_version = quast.version
    Int? assembly_length = quast.genome_length
    Int? number_contigs = quast.number_contigs
    Int? n50_value = quast.n50_value
    Float? quast_gc_percent = quast.gc_percent
    # Assembly QC - nanoplot outputs
    Float? est_coverage_raw = nanoplot_raw.est_coverage
    Float? est_coverage_clean = nanoplot_clean.est_coverage
    # Assembly QC - busco outputs
    String? busco_version = busco.busco_version
    String? busco_docker = busco.busco_docker
    String? busco_database = busco.busco_database
    String? busco_results = busco.busco_results
    File? busco_report = busco.busco_report
    # Gambit outputs
    File? gambit_report_file = gambit.gambit_report_file
    File? gambit_closest_genomes_file = gambit.gambit_closest_genomes_file
    String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String? gambit_next_taxon = gambit.gambit_next_taxon
    String? gambit_next_taxon_rank = gambit.gambit_next_taxon_rank
    String? gambit_version = gambit.gambit_version
    String? gambit_db_version = gambit.gambit_db_version
    String? merlin_tag = gambit.merlin_tag
    String? gambit_docker = gambit.gambit_docker
    # C. auris specific outputs for cladetyper
    String? cladetyper_clade = merlin_magic.clade_type
    String? cladetyper_version = merlin_magic.cladetyper_version
    String? cladetyper_docker_image = merlin_magic.cladetyper_docker_image
    String? cladetype_annotated_ref = merlin_magic.cladetype_annotated_ref
    # AMR Search outputs
    File? amr_search_results = merlin_magic.amr_search_results
    File? amr_search_csv = merlin_magic.amr_results_csv
    File? amr_search_results_pdf = merlin_magic.amr_results_pdf
    String? amr_search_docker = merlin_magic.amr_search_docker
    String? amr_search_version = merlin_magic.amr_search_version  
    # Snippy variants outputs
    String? theiaeuk_snippy_variants_version = merlin_magic.snippy_variants_version
    String? theiaeuk_snippy_variants_query = merlin_magic.snippy_variants_query
    String? theiaeuk_snippy_variants_query_check = merlin_magic.snippy_variants_query_check
    String? theiaeuk_snippy_variants_hits = merlin_magic.snippy_variants_hits
    String? theiaeuk_snippy_variants_gene_query_results = merlin_magic.snippy_variants_gene_query_results
    String? theiaeuk_snippy_variants_outdir_tarball = merlin_magic.snippy_variants_outdir_tarball
    String? theiaeuk_snippy_variants_results = merlin_magic.snippy_variants_results
    String? theiaeuk_snippy_variants_bam = merlin_magic.snippy_variants_bam
    String? theiaeuk_snippy_variants_bai = merlin_magic.snippy_variants_bai
    String? theiaeuk_snippy_variants_summary = merlin_magic.snippy_variants_summary
    String? theiaeuk_snippy_variants_num_reads_aligned = merlin_magic.snippy_variants_num_reads_aligned
    String? theiaeuk_snippy_variants_coverage_tsv = merlin_magic.snippy_variants_coverage_tsv
    String? theiaeuk_snippy_variants_num_variants = merlin_magic.snippy_variants_num_variants
    String? theiaeuk_snippy_variants_percent_ref_coverage = merlin_magic.snippy_variants_percent_ref_coverage
  }
}