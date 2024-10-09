version 1.0

import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc
import "../../tasks/assembly/task_dragonflye.wdl" as dragonflye
import "../../tasks/taxon_id/task_gambit.wdl" as gambit
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/quality_control/basic_statistics/task_nanoplot.wdl" as nanoplot_task
import "../../tasks/quality_control/basic_statistics/task_quast.wdl" as quast_task
import "../../tasks/quality_control/advanced_metrics/task_busco.wdl" as busco_task

workflow theiaeuk_ont {
  input {
    File read1
    String samplename
    Int genome_length = 50000000 
    String workflow_series = "theiaeuk"
    String? assembler
    String? assembler_options
    Int dragonflye_cpu = 8
    Int dragonflye_memory = 32
    Int dragonflye_disk_size = 100
    String medaka_model = "r941_min_hac_g507"
    File gambit_db_genomes = "gs://theiagen-public-files-rp/terra/theiaeuk-files/gambit/221130-theiagen-fungal-v0.2.db"
    File gambit_db_signatures = "gs://theiagen-public-files-rp/terra/theiaeuk-files/gambit/221130-theiagen-fungal-v0.2.h5"
  }

  call versioning.version_capture {
    input:
  }

  call read_qc.read_QC_trim_ont as read_qc {
    input:
      read1 = read1,
      samplename = samplename,
      genome_length = genome_length,
      workflow_series = workflow_series
  }

  call dragonflye.dragonflye {
    input:
      read1 = read_qc.read1_clean,
      samplename = samplename,
      assembler = assembler,
      assembler_options = assembler_options,
      genome_length = genome_length,
      medaka_model = medaka_model,
      cpu = dragonflye_cpu,
      memory = dragonflye_memory,
      disk_size = dragonflye_disk_size
  }
  #call quast on the assembly
  call quast_task.quast {
    input:
      assembly = dragonflye.assembly_fasta,
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
      read1 = read_qc.read1_clean,
      samplename = samplename,
      est_genome_length = select_first([genome_length, quast.genome_length])
  }
  # busco on the assembly
  call busco_task.busco {
    input:
      assembly = dragonflye.assembly_fasta,
      samplename = samplename
  }

  call gambit.gambit {
    input:
      assembly = dragonflye.assembly_fasta,
      samplename = samplename,
      gambit_db_genomes = gambit_db_genomes,
      gambit_db_signatures = gambit_db_signatures
  }

   call merlin_magic_workflow.merlin_magic {
    input:
      samplename = samplename,
      merlin_tag = gambit.merlin_tag,
      assembly = dragonflye.assembly_fasta,
      assembly_only = true,
      ont_data = true,
      theiaeuk = true
  }

  output {
    # Version Capture
    String theiaeuk_ont_version = version_capture.phb_version
    String theiaeuk_ont_analysis_date = version_capture.date
    # Read QC outputs
    File read1_clean = read_qc.read1_clean
    String? nanoq_version = read_qc.nanoq_version
    Int est_genome_length = read_qc.est_genome_length
    # Assembly outputs
    File assembly_fasta = dragonflye.assembly_fasta
    File contigs_gfa = dragonflye.contigs_gfa
    String dragonflye_version = dragonflye.dragonflye_version
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
    File gambit_report_file = gambit.gambit_report_file
    File gambit_closest_genomes_file = gambit.gambit_closest_genomes_file
    String gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String gambit_next_taxon = gambit.gambit_next_taxon
    String gambit_next_taxon_rank = gambit.gambit_next_taxon_rank
    String gambit_version = gambit.gambit_version
    String gambit_db_version = gambit.gambit_db_version
    String merlin_tag = gambit.merlin_tag
    String gambit_docker = gambit.gambit_docker
    # C. auris specific outputs for cladetyper
    String? clade_type = merlin_magic.clade_type
    String? cladetyper_analysis_date = merlin_magic.cladetyper_analysis_date
    String? cladetyper_version = merlin_magic.cladetyper_version
    String? cladetyper_docker_image = merlin_magic.cladetyper_docker_image
    String? cladetype_annotated_ref = merlin_magic.cladetype_annotated_ref
    # Snippy variants outputs
    String? snippy_variants_version = merlin_magic.snippy_variants_version
    String? snippy_variants_query = merlin_magic.snippy_variants_query
    String? snippy_variants_query_check = merlin_magic.snippy_variants_query_check
    String? snippy_variants_hits = merlin_magic.snippy_variants_hits
    String? snippy_variants_gene_query_results = merlin_magic.snippy_variants_gene_query_results
    String? snippy_variants_outdir_tarball = merlin_magic.snippy_variants_outdir_tarball
    String? snippy_variants_results = merlin_magic.snippy_variants_results
    String? snippy_variants_bam = merlin_magic.snippy_variants_bam
    String? snippy_variants_bai = merlin_magic.snippy_variants_bai
    String? snippy_variants_summary = merlin_magic.snippy_variants_summary
    String? snippy_variants_num_reads_aligned = merlin_magic.snippy_variants_num_reads_aligned
    String? snippy_variants_coverage_tsv = merlin_magic.snippy_variants_coverage_tsv
    String? snippy_variants_num_variants = merlin_magic.snippy_variants_num_variants
    String? snippy_variants_percent_ref_coverage = merlin_magic.snippy_variants_percent_ref_coverage
  }
}