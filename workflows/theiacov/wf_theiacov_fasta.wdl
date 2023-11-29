
version 1.0

import "../../tasks/quality_control/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/gene_typing/task_abricate.wdl" as abricate
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/species_typing/task_pangolin.wdl" as pangolin
import "../../tasks/quality_control/task_qc_check_phb.wdl" as qc_check
import "../../tasks/utilities/task_organism_defaults.wdl" as defaults
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
  if (organism == "sars-cov-2") {
    call defaults.set_organism_defaults_sc2 {
      input:
        reference_genome = reference_genome, 
        nextclade_ref = nextclade_dataset_reference,
        nextclade_ds_tag = nextclade_dataset_tag,
        nextclade_ds_name = nextclade_dataset_name,
        genome_len = genome_length,
        vadr_max_length = maxlen, 
        vadr_options = vadr_opts
    }
  }
  if (organism == "MPXV") {
    call defaults.set_organism_defaults_mpox {
      input:
        reference_genome = reference_genome,
        nextclade_ref = nextclade_dataset_reference,
        nextclade_ds_tag = nextclade_dataset_tag,
        nextclade_ds_name = nextclade_dataset_name,
        genome_len = genome_length,
        vadr_max_length = maxlen, 
        vadr_options = vadr_opts
    }
  }
  if (organism == "WNV") {
    call defaults.set_organism_defaults_wnv {
      input:
        reference_genome = reference_genome,
        genome_len = genome_length,
        vadr_max_length = maxlen, 
        vadr_options = vadr_opts
    }
  }
  if (organism == "flu") {
    if (!defined(flu_subtype)) {
      call abricate.abricate_flu {
        input:
          assembly = assembly_fasta,
          samplename = samplename
      }
      String? abricate_subtype = abricate_flu.abricate_flu_subtype
    }
    call defaults.set_organism_defaults_flu {
      input:
        flu_segment = flu_segment,
        flu_subtype = select_first([flu_subtype, abricate_subtype]),
        genome_len = genome_length
    }
  }
  if (organism == "rsv_a") {
    call defaults.set_organism_defaults_rsv_a {
      input:
        reference_genome = reference_genome,
        nextclade_ds_tag = nextclade_dataset_tag,
        nextclade_ds_name = nextclade_dataset_name,
        nextclade_ref = nextclade_dataset_reference,
        genome_len = genome_length,
        vadr_max_length = maxlen, 
        vadr_options = vadr_opts
    }
  } 
  if (organism == "rsv_b") {
    call defaults.set_organism_defaults_rsv_b {
      input:
        reference_genome = reference_genome,
        nextclade_ds_tag = nextclade_dataset_tag,
        nextclade_ds_name = nextclade_dataset_name,
        nextclade_ref = nextclade_dataset_reference,
        genome_len = genome_length,
        vadr_max_length = maxlen, 
        vadr_options = vadr_opts
    }
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = select_first([reference_genome, set_organism_defaults_sc2.reference, set_organism_defaults_mpox.reference, set_organism_defaults_wnv.reference, set_organism_defaults_flu.reference, set_organism_defaults_rsv_a.reference, set_organism_defaults_rsv_b.reference]),
      genome_length = select_first([set_organism_defaults_sc2.genome_length, set_organism_defaults_mpox.genome_length, set_organism_defaults_wnv.genome_length, set_organism_defaults_flu.genome_length, set_organism_defaults_rsv_a.genome_length, set_organism_defaults_rsv_b.genome_length]),
  }
  if (organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta
    }
  }
  if (organism == "sars-cov-2" || organism == "MPXV" || organism == "rsv_a" || organism == "rsv_b" || organism == "flu") {
    if (select_first([set_organism_defaults_flu.nextclade_dataset_tag, ""]) != "NA") {
      call nextclade_task.nextclade {
        input:
          genome_fasta = assembly_fasta,
          dataset_name = select_first([nextclade_dataset_name, set_organism_defaults_sc2.nextclade_dataset_name, set_organism_defaults_mpox.nextclade_dataset_name, set_organism_defaults_rsv_a.nextclade_dataset_name, set_organism_defaults_rsv_b.nextclade_dataset_name, set_organism_defaults_flu.nextclade_dataset_name]),
          dataset_reference = select_first([nextclade_dataset_reference, set_organism_defaults_sc2.nextclade_reference, set_organism_defaults_mpox.nextclade_reference, set_organism_defaults_rsv_a.nextclade_reference, set_organism_defaults_rsv_b.nextclade_reference, set_organism_defaults_flu.nextclade_reference]),
          dataset_tag = select_first([nextclade_dataset_tag, set_organism_defaults_sc2.nextclade_dataset_tag, set_organism_defaults_mpox.nextclade_dataset_tag, set_organism_defaults_rsv_a.nextclade_dataset_tag, set_organism_defaults_rsv_b.nextclade_dataset_tag, set_organism_defaults_flu.nextclade_dataset_tag])
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
        maxlen = select_first([maxlen, set_organism_defaults_sc2.vadr_maxlen, set_organism_defaults_mpox.vadr_maxlen, set_organism_defaults_rsv_a.vadr_maxlen, set_organism_defaults_rsv_b.vadr_maxlen, set_organism_defaults_wnv.vadr_maxlen]),
        vadr_opts = select_first([vadr_opts, set_organism_defaults_sc2.vadr_opts, set_organism_defaults_mpox.vadr_opts, set_organism_defaults_rsv_a.vadr_opts, set_organism_defaults_rsv_b.vadr_opts, set_organism_defaults_wnv.vadr_opts]),
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
    String nextclade_ds_tag =  select_first([nextclade_dataset_tag, set_organism_defaults_sc2.nextclade_dataset_tag, set_organism_defaults_mpox.nextclade_dataset_tag, set_organism_defaults_rsv_a.nextclade_dataset_tag, set_organism_defaults_rsv_b.nextclade_dataset_tag, set_organism_defaults_flu.nextclade_dataset_tag, "NA"])
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
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