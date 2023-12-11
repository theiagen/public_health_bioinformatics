
version 1.0

import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/task_versioning.wdl" as versioning

workflow theiacov_fasta {
  meta {
    description: "Reference-based consensus calling for viral amplicon ont sequencing data generated on the Clear Labs platform."
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2"
    # sequencing values
    String seq_method
    String input_assembly_method
    # nextclade inputs
    String nextclade_dataset_reference = "MN908947"
    String nextclade_dataset_tag = "2023-09-21T12:00:00Z"
    String? nextclade_dataset_name
    # qc check parameters
    File? qc_check_table
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta
  }
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta
    }
  }
  if (organism == "MPXV") {
    # MPXV specific tasks
  }
  if (organism == "WNV") {
    # WNV specific tasks (none yet, just adding as placeholder for future)
  }
  if (organism == "MPXV" || organism == "sars-cov-2"){
    # tasks specific to either MPXV or sars-cov-2 
    call nextclade_task.nextclade {
      input:
      genome_fasta = assembly_fasta,
      dataset_name = select_first([nextclade_dataset_name, organism]),
      dataset_reference = nextclade_dataset_reference,
      dataset_tag = nextclade_dataset_tag
    }
    call nextclade_task.nextclade_output_parser {
      input:
      nextclade_tsv = nextclade.nextclade_tsv,
      organism = organism
    }
  }
  if (organism == "MPXV" || organism == "sars-cov-2" || organism == "WNV"){ 
    # tasks specific to MPXV, sars-cov-2, and WNV
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG
    }
  }
  if(defined(qc_check_table)) {
    call qc_check.qc_check_phb as qc_check_task {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism,
        number_N = consensus_qc.number_N,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        number_Degenerate =  consensus_qc.number_Degenerate,
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
    String nextclade_ds_tag = nextclade_dataset_tag
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
  }
}