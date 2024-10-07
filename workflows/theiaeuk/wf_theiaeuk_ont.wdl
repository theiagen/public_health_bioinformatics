version 1.0

import "../utilities/wf_read_QC_trim_ont.wdl" as read_qc
import "../../tasks/assembly/task_dragonflye.wdl" as dragonflye
import "../../tasks/taxon_id/task_gambit.wdl" as gambit
import "../utilities/wf_merlin_magic.wdl" as merlin_magic_workflow
import "../../tasks/task_versioning.wdl" as versioning

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
    Boolean assembly_only = true 
    Boolean ont_data = true
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
      read1 = read_qc.read1_clean,
      assembly_only = assembly_only,
      ont_data = ont_data,
      theiaeuk = true
  }

  call versioning.version_capture {
    input:
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