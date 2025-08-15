version 1.0

import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check_phb_task
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_flu_track.wdl" as run_flu_track

workflow theiacov_fasta {
  meta {
    description: "Assessment of the quality of a consensus assembly fasta file for sars-cov-2, MPXV, WNV, flu, or RSV"
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b
    # flu options
    String? flu_segment # options: HA or NA
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1" "H5N1"
    # optional reference information
    File? reference_genome
    Int? genome_length
    # nextclade inputs (default SC2)
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # sequencing values
    String seq_method
    String input_assembly_method
    # qc check parameters
    File? qc_check_table
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_opts
    File? vadr_model_file
    Int? vadr_memory
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      flu_segment = flu_segment,
      flu_subtype = flu_subtype,
      reference_genome = reference_genome,
      genome_length_input = genome_length,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,
      vadr_max_length = vadr_max_length,
      vadr_skip_length = vadr_skip_length,
      vadr_options = vadr_opts,
      vadr_model = vadr_model_file,
      vadr_mem = vadr_memory
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = organism_parameters.reference,
      genome_length = organism_parameters.genome_length
  }
  # vadr task
if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles") {
  call vadr_task.vadr {
    input:
      genome_fasta = assembly_fasta,
      assembly_length_unambiguous = consensus_qc.number_ATCG,
      max_length = organism_parameters.vadr_maxlength,
      vadr_opts = organism_parameters.vadr_opts,
      vadr_model_file = organism_parameters.vadr_model_file,
      skip_length = organism_parameters.vadr_skiplength,
      memory = organism_parameters.vadr_memory
    }
  }
  if (organism == "flu") {
    call run_flu_track.flu_track {
      input:
        assembly_fasta = assembly_fasta,
        samplename = samplename,
        standardized_organism = organism,
        seq_method = seq_method,
        flu_subtype = flu_subtype,
        vadr_outputs_tgz = vadr.outputs_tgz,
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta,
        docker = organism_parameters.pangolin_docker
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "measles") {
    if (organism_parameters.nextclade_dataset_tag != "NA") {
      call nextclade_task.nextclade_v3 {
        input:
          genome_fasta = assembly_fasta,
          dataset_name = organism_parameters.nextclade_dataset_name,
          dataset_tag = organism_parameters.nextclade_dataset_tag
      }
      call nextclade_task.nextclade_output_parser {
        input:
          nextclade_tsv = nextclade_v3.nextclade_tsv,
          organism = organism_parameters.standardized_organism
      }
    }
  }
  # QC check task
  if (defined(qc_check_table)) {
    call qc_check_phb_task.qc_check_phb as qc_check_task {
      input:
        qc_check_table = qc_check_table,
        expected_taxon = organism_parameters.standardized_organism,
        number_N = consensus_qc.number_N,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        number_Degenerate = consensus_qc.number_Degenerate,
        percent_reference_coverage =  consensus_qc.percent_reference_coverage,
        vadr_num_alerts = vadr.num_alerts
    }
  }
  call versioning.version_capture {
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
    File? nextclade_json = nextclade_v3.nextclade_json
    File? auspice_json = nextclade_v3.auspice_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    String? nextclade_version = nextclade_v3.nextclade_version
    String? nextclade_docker = nextclade_v3.nextclade_docker
    String nextclade_ds_tag =  organism_parameters.nextclade_dataset_tag
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # VADR Annotation QC
    File? vadr_alerts_list = vadr.alerts_list
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
    String? vadr_num_alerts = vadr.num_alerts
    # VADR Annotation QC for flu
    File? vadr_flu_segment_concatenated_fasta = flu_track.flu_assembly_fasta_concatenated
    File? vadr_flu_ha_segment_fasta = flu_track.flu_ha_segment_fasta
    File? vadr_flu_na_segment_fasta = flu_track.flu_na_segment_fasta
    File? vadr_flu_pa_segment_fasta = flu_track.flu_pa_segment_fasta
    File? vadr_flu_pb1_segment_fasta = flu_track.flu_pb1_segment_fasta
    File? vadr_flu_pb2_segment_fasta = flu_track.flu_pb2_segment_fasta
    File? vadr_flu_mp_segment_fasta = flu_track.flu_mp_segment_fasta
    File? vadr_flu_np_segment_fasta = flu_track.flu_np_segment_fasta
    File? vadr_flu_ns_segment_fasta = flu_track.flu_ns_segment_fasta
    # QC_Check Results
    String? qc_check = qc_check_task.qc_check
    File? qc_standard = qc_check_task.qc_standard
    # Flu Outputs
    String? abricate_flu_type = flu_track.abricate_flu_type
    String? abricate_flu_subtype =  flu_track.abricate_flu_subtype
    File? abricate_flu_results = flu_track.abricate_flu_results
    String? abricate_flu_database =  flu_track.abricate_flu_database
    String? abricate_flu_version = flu_track.abricate_flu_version
    # GenoFLU outputs
    String? genoflu_version = flu_track.genoflu_version
    String? genoflu_genotype = flu_track.genoflu_genotype
    String? genoflu_all_segments = flu_track.genoflu_all_segments
    File? genoflu_output_tsv = flu_track.genoflu_output_tsv
    # Nextclade outputs for flu H5N1
    File? nextclade_json_flu_h5n1 = flu_track.nextclade_json_flu_h5n1
    File? auspice_json_flu_h5n1 = flu_track.auspice_json_flu_h5n1
    File? nextclade_tsv_flu_h5n1 = flu_track.nextclade_tsv_flu_h5n1
    String? nextclade_aa_subs_flu_h5n1 = flu_track.nextclade_aa_subs_flu_h5n1
    String? nextclade_aa_dels_flu_h5n1 = flu_track.nextclade_aa_dels_flu_h5n1
    String? nextclade_clade_flu_h5n1 = flu_track.nextclade_clade_flu_h5n1
    String? nextclade_qc_flu_h5n1 = flu_track.nextclade_qc_flu_h5n1
    # Nextclade outputs for flu HA
    File? nextclade_json_flu_ha = flu_track.nextclade_json_flu_ha
    File? auspice_json_flu_ha = flu_track.auspice_json_flu_ha
    File? nextclade_tsv_flu_ha = flu_track.nextclade_tsv_flu_ha
    String? nextclade_ds_tag_flu_ha = flu_track.nextclade_ds_tag_flu_ha
    String? nextclade_aa_subs_flu_ha = flu_track.nextclade_aa_subs_flu_ha
    String? nextclade_aa_dels_flu_ha = flu_track.nextclade_aa_dels_flu_ha
    String? nextclade_clade_flu_ha = flu_track.nextclade_clade_flu_ha
    String? nextclade_qc_flu_ha = flu_track.nextclade_qc_flu_ha
    # Nextclade outputs for flu NA
    File? nextclade_json_flu_na = flu_track.nextclade_json_flu_na
    File? auspice_json_flu_na = flu_track.auspice_json_flu_na
    File? nextclade_tsv_flu_na = flu_track.nextclade_tsv_flu_na
    String? nextclade_ds_tag_flu_na = flu_track.nextclade_ds_tag_flu_na
    String? nextclade_aa_subs_flu_na = flu_track.nextclade_aa_subs_flu_na
    String? nextclade_aa_dels_flu_na = flu_track.nextclade_aa_dels_flu_na
    String? nextclade_clade_flu_na = flu_track.nextclade_clade_flu_na
    String? nextclade_qc_flu_na = flu_track.nextclade_qc_flu_na
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = flu_track.flu_A_315675_resistance
    String? flu_amantadine_resistance = flu_track.flu_amantadine_resistance
    String? flu_compound_367_resistance = flu_track.flu_compound_367_resistance
    String? flu_favipiravir_resistance = flu_track.flu_favipiravir_resistance
    String? flu_fludase_resistance = flu_track.flu_fludase_resistance
    String? flu_L_742_001_resistance = flu_track.flu_L_742_001_resistance
    String? flu_laninamivir_resistance = flu_track.flu_laninamivir_resistance
    String? flu_peramivir_resistance = flu_track.flu_peramivir_resistance
    String? flu_pimodivir_resistance = flu_track.flu_pimodivir_resistance
    String? flu_rimantadine_resistance = flu_track.flu_rimantadine_resistance
    String? flu_oseltamivir_resistance = flu_track.flu_oseltamivir_resistance
    String? flu_xofluza_resistance = flu_track.flu_xofluza_resistance
    String? flu_zanamivir_resistance = flu_track.flu_zanamivir_resistance
  }
}