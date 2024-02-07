
version 1.0

import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../workflows/utilities/wf_organism_parameters.wdl" as defaults
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta {
  meta {
    description: "Assessment of the quality of a consensus assembly fasta file for sars-cov-2, MPXV, WNV, flu, or RSV"
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b
    # flu options
    String flu_segment = "HA" # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1"
    # optional reference information
    File? reference_genome
    Int? genome_length
    # nextclade inputs (default SC2)
    String? nextclade_dataset_reference
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # sequencing values
    String seq_method
    String input_assembly_method
    # qc check parameters
    File? qc_check_table
    # vadr parameters
    Int? maxlen
    String? vadr_opts
  }
  # only run abricate if user sets organism = "flu" AND if flu_subtype is unknown/not set by user
  if (!defined(flu_subtype) && organism == "flu") {
    call abricate.abricate_flu {
      input:
        assembly = assembly_fasta,
        samplename = samplename
    }
  String abricate_subtype = abricate_flu.abricate_flu_subtype
  }
  call defaults.organism_parameters {
    input:
      organism = organism,
      flu_segment = flu_segment,
      flu_subtype = select_first([flu_subtype, abricate_subtype, "N/A"]),
      reference_genome = reference_genome,
      genome_length = genome_length,
      nextclade_ds_reference = nextclade_dataset_reference,
      nextclade_ds_tag = nextclade_dataset_tag,
      nextclade_ds_name = nextclade_dataset_name,
      vadr_max_length = maxlen,
      vadr_options = vadr_opts
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = organism_parameters.reference,
      genome_length = organism_parameters.genome_len
  }
  if (organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta
    }
  }
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (organism_parameters.nextclade_dataset_tag != "NA") {
      call nextclade_task.nextclade {
        input:
          genome_fasta = assembly_fasta,
          dataset_name = organism_parameters.nextclade_dataset_name,
          dataset_reference = organism_parameters.nextclade_dataset_reference,
          dataset_tag = organism_parameters.nextclade_dataset_tag
      }
    }
  }
  # nextclade parser task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (defined(nextclade.nextclade_tsv)) {
      call nextclade_task.nextclade_output_parser {
        input:
          nextclade_tsv = select_first([nextclade.nextclade_tsv]),
          organism = organism
      }
    }
  }
  # vadr task
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "WNV") {
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        maxlen = organism_parameters.vadr_maxlen,
        vadr_opts = organism_parameters.vadr_opts
    }
  }
  # QC check task
  if(defined(qc_check_table)) {
    call qc_check.qc_check_phb {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism,
        number_N = consensus_qc.number_N,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        number_Degenerate = consensus_qc.number_Degenerate,
        percent_reference_coverage =  consensus_qc.percent_reference_coverage,
        vadr_num_alerts = vadr.num_alerts
    }
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String theiacov_fasta_version = version_capture.phb_version
    String theiacov_fasta_analysis_date = version_capture.date
    # Read & Assembly Metadata
    String seq_platform = seq_method
    String assembly_method = input_assembly_method
    # Assembly QC - consensus assembly summary statistics
    Int number_N = consensus_qc.number_N
    Int assembly_length_unambiguous = consensus_qc.number_ATCG
    Int number_Degenerate = consensus_qc.number_Degenerate
    Int number_Total = consensus_qc.number_Total
    Float percent_reference_coverage = consensus_qc.percent_reference_coverage
    # Pangolin outputs
    String? pango_lineage = pangolin4.pangolin_lineage
    String? pango_lineage_expanded = pangolin4.pangolin_lineage_expanded
    String? pangolin_conflicts = pangolin4.pangolin_conflicts
    String? pangolin_notes = pangolin4.pangolin_notes
    String? pangolin_assignment_version = pangolin4.pangolin_assignment_version
    File? pango_lineage_report = pangolin4.pango_lineage_report
    String? pangolin_docker = pangolin4.pangolin_docker
    String? pangolin_versions = pangolin4.pangolin_versions
    # Nextclade outputs
    File? nextclade_json = nextclade.nextclade_json
    File? auspice_json = nextclade.auspice_json
    File? nextclade_tsv = nextclade.nextclade_tsv
    String? nextclade_version = nextclade.nextclade_version
    String? nextclade_docker = nextclade.nextclade_docker
    String nextclade_ds_tag =  organism_parameters.nextclade_dataset_tag
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    String? vadr_num_alerts = vadr.num_alerts
    # QC_Check Results
    String? qc_check = qc_check_phb.qc_check
    File? qc_standard = qc_check_phb.qc_standard
    # Flu Outputs
    String? abricate_flu_type = abricate_flu.abricate_flu_type
    String? abricate_flu_subtype =  abricate_flu.abricate_flu_subtype
    File? abricate_flu_results = abricate_flu.abricate_flu_results
    String? abricate_flu_database =  abricate_flu.abricate_flu_database
    String? abricate_flu_version = abricate_flu.abricate_flu_version
  }
}