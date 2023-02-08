version 1.0

import "wf_read_QC_trim_pe_theiaeuk.wdl" as read_qc
import "wf_merlin_magic_euk.wdl" as merlin_magic
import "../tasks/assembly/task_shovill.wdl" as shovill
import "../tasks/quality_control/task_quast.wdl" as quast
import "../tasks/quality_control/task_cg_pipeline.wdl" as cg_pipeline
import "../tasks/quality_control/task_screen.wdl" as screen
import "../tasks/taxon_id/task_gambit.wdl" as gambit
import "../tasks/species_typing/task_ts_mlst.wdl" as ts_mlst
import "../tasks/task_versioning.wdl" as versioning

workflow theiaeuk_illumina_pe {
  meta {
    description: "De-novo genome assembly, taxonomic ID, and QC of paired-end eukaryotic NGS data"
  }
  input {
    String samplename
    String seq_method = "ILLUMINA"
    File read1_raw
    File read2_raw
    Int min_reads = 30000
    #Edited default values
    Int min_basepairs = 90000000
    Int min_genome_size = 9000000
    Int max_genome_size = 178000000
    Int min_coverage = 10
    Int min_proportion = 50
    Boolean skip_screen = false 
  }
  call versioning.version_capture{
    input:
  }
  call screen.check_reads as raw_check_reads {
    input:
      read1 = read1_raw,
      read2 = read2_raw,
      min_reads = min_reads,
      min_basepairs = min_basepairs,
      min_genome_size = min_genome_size,
      max_genome_size = max_genome_size,
      min_coverage = min_coverage,
      min_proportion = min_proportion,
      skip_screen = skip_screen
  }
  if (raw_check_reads.read_screen=="PASS") {
    call read_qc.read_QC_trim_pe as read_QC_trim {
      input:
        samplename = samplename,
        read1_raw = read1_raw,
        read2_raw = read2_raw
    }
    call screen.check_reads as clean_check_reads {
      input:
        read1 = read_QC_trim.read1_clean,
        read2 = read_QC_trim.read2_clean,
        min_reads = min_reads,
        min_basepairs = min_basepairs,
        min_genome_size = min_genome_size,
        max_genome_size = max_genome_size,
        min_coverage = min_coverage,
        min_proportion = min_proportion,
        skip_screen = skip_screen
    }
    if (clean_check_reads.read_screen=="PASS") {
      call shovill.shovill_pe {
        input:
          samplename = samplename,
          read1_cleaned = read_QC_trim.read1_clean,
          read2_cleaned = read_QC_trim.read2_clean
      }
      call quast.quast {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call cg_pipeline.cg_pipeline {
        input:
          read1 = read1_raw,
          read2 = read2_raw,
          samplename = samplename,
          genome_length = clean_check_reads.est_genome_length
      }
      call gambit.gambit_euk as gambit {
        input:
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call ts_mlst.ts_mlst {
        input: 
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename
      }
      call merlin_magic.merlin_magic {
        input:
          merlin_tag = gambit.merlin_tag,
          assembly = shovill_pe.assembly_fasta,
          samplename = samplename,
          read1 = read_QC_trim.read1_clean,
          read2 = read_QC_trim.read2_clean
      }
    }
  }
  output {
    # Version Captures
    String theiaeuk_illumina_pe_version = version_capture.phbg_version
    String theiaeuk_illumina_pe_analysis_date = version_capture.date
    # Read Metadata
    String seq_platform = seq_method
    # Sample Screening
    String raw_read_screen = raw_check_reads.read_screen
    String? clean_read_screen = clean_check_reads.read_screen
    # Read QC
    Int? num_reads_raw1 = read_QC_trim.fastq_scan_raw1
    Int? num_reads_raw2 = read_QC_trim.fastq_scan_raw2
    String? num_reads_raw_pairs = read_QC_trim.fastq_scan_raw_pairs
    String? fastq_scan_version = read_QC_trim.fastq_scan_version
    Int? num_reads_clean1 = read_QC_trim.fastq_scan_clean1
    Int? num_reads_clean2 = read_QC_trim.fastq_scan_clean2
    String? num_reads_clean_pairs = read_QC_trim.fastq_scan_clean_pairs
    String? trimmomatic_version = read_QC_trim.trimmomatic_version
    String? bbduk_docker = read_QC_trim.bbduk_docker
    Float? r1_mean_q = cg_pipeline.r1_mean_q
    Float? r2_mean_q = cg_pipeline.r2_mean_q
    File? read1_clean = read_QC_trim.read1_clean
    File? read2_clean = read_QC_trim.read2_clean
    #Assembly and Assembly QC
    File? assembly_fasta = shovill_pe.assembly_fasta
    File? contigs_gfa = shovill_pe.contigs_gfa
    File? contigs_fastg = shovill_pe.contigs_fastg
    File? contigs_lastgraph = shovill_pe.contigs_lastgraph
    String? shovill_pe_version = shovill_pe.shovill_version
    File? quast_report = quast.quast_report
    String? quast_version = quast.version
    Int? genome_length = quast.genome_length
    Int? number_contigs = quast.number_contigs
    Int? n50_value = quast.n50_value
    File? cg_pipeline_report = cg_pipeline.cg_pipeline_report
    String? cg_pipeline_docker = cg_pipeline.cg_pipeline_docker
    Float? est_coverage = cg_pipeline.est_coverage
    # Taxon ID
    File? gambit_report = gambit.gambit_report_file
    File? gambit_closest_genomes = gambit.gambit_closest_genomes_file
    String? gambit_predicted_taxon = gambit.gambit_predicted_taxon
    String? gambit_predicted_taxon_rank = gambit.gambit_predicted_taxon_rank
    String? gambit_version = gambit.gambit_version
    String? gambit_db_version = gambit.gambit_db_version
    String? gambit_docker = gambit.gambit_docker
    # MLST Typing
    File? ts_mlst_results = ts_mlst.ts_mlst_results
    String? ts_mlst_predicted_st = ts_mlst.ts_mlst_predicted_st
    String? ts_mlst_version = ts_mlst.ts_mlst_version
    String? ts_mlst_pubmlst_scheme = ts_mlst.ts_mlst_pubmlst_scheme
    # Cladetyper Outputs
    String? clade_type = merlin_magic.clade_type
    String? cladetyper_analysis_date = merlin_magic.cladetyper_analysis_date
    String? cladetyper_version = merlin_magic.cladetyper_version
    String? cladetyper_docker_image = merlin_magic.cladetyper_docker_image
    String? cladetype_annotated_ref = merlin_magic.cladetype_annotated_ref
    # Snippy Outputs
    String? theiaeuk_snippy_variants_version = merlin_magic.snippy_variants_version
    String? theiaeuk_snippy_variants_query = merlin_magic.snippy_variants_query
    String? theiaeuk_snippy_variants_hits = merlin_magic.snippy_variants_hits
    File? theiaeuk_snippy_variants_gene_query_results = merlin_magic.snippy_variants_gene_query_results
    # Array[File]? snippy_outputs = merlin_magic.snippy_outputs
    File? theiaeuk_snippy_variants_results = merlin_magic.snippy_variants_results
    File? theiaeuk_snippy_variants_bam = merlin_magic.snippy_variants_bam
    File? theiaeuk_snippy_variants_bai = merlin_magic.snippy_variants_bai
    File? theiaeuk_snippy_variants_summary = merlin_magic.snippy_variants_summary
  }
}