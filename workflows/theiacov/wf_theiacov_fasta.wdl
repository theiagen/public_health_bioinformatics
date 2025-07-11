
version 1.0

import "../../tasks/gene_typing/drug_resistance/task_abricate.wdl" as abricate
import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_consensus_qc.wdl" as consensus_qc_task
import "../../tasks/quality_control/comparisons/task_qc_check_phb.wdl" as qc_check
import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/species_typing/orthomyxoviridae/task_genoflu.wdl" as genoflu_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_influenza_antiviral_substitutions.wdl" as flu_antiviral
import "../../tasks/species_typing/influenza/task_extract_flu_segments.wdl" as extract_flu_segments_task

workflow theiacov_fasta {
  meta {
    description: "Assessment of the quality of a consensus assembly fasta file for sars-cov-2, MPXV, WNV, flu, or RSV"
  }
  input {
    String samplename
    File assembly_fasta
    String organism = "sars-cov-2" # options: "sars-cov-2" "MPXV" "WNV" "flu" "rsv_a" "rsv_b
    # flu options
    String? flu_segment # options: HA or NA. only required if input assembly is a singular flu segment.
    String? flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1" "H5N1"
    # optional reference information
    File? reference_genome
    Int? genome_length
    # Abricate inputs
    Int? abricate_flu_min_percent_identity
    Int? abricate_flu_min_percent_coverage
    String? abricate_flu_docker
    Int? abricate_flu_memory
    Int? abricate_flu_cpu
    Int? abricate_flu_disk_size
    # flu antiviral substitutions subworkflow inputs
    File? flu_h1_ha_ref
    File? flu_h3_ha_ref
    File? flu_n1_na_ref
    File? flu_n2_na_ref
    File? flu_pa_ref
    File? flu_pb1_ref
    File? flu_pb2_ref
    File? flu_h1n1_m2_ref
    File? flu_h3n2_m2_ref
    String? antiviral_aa_subs
    # nextclade inputs (default SC2)
    String? nextclade_docker_image
    Int? nextclade_cpu
    Int? nextclade_memory
    Int? nextclade_disk_size
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    File? nextclade_custom_input_dataset
    # nextclade output parser inputs
    String? nextclade_output_parser_docker
    Int? nextclade_output_parser_cpu
    Int? nextclade_output_parser_memory
    Int? nextclade_output_parser_disk_size
    # sequencing values
    String seq_method
    String input_assembly_method
    # qc check parameters
    File? qc_check_table
    # GenoFLU inputs
    Float? genoflu_min_percent_identity
    File? genoflu_cross_reference
    Int? genoflu_cpu
    Int? genoflu_disk_size
    String? genoflu_docker
    Int? genoflu_memory
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_opts
    Int? vadr_memory
  }
  # only run abricate if user sets organism = "flu" AND if flu_subtype is unknown/not set by user
  if (organism == "flu") {
    call abricate.abricate_flu {
      input:
        assembly = assembly_fasta,
        samplename = samplename,
        min_percent_identity = abricate_flu_min_percent_identity,
        min_percent_coverage = abricate_flu_min_percent_coverage,
        cpu = abricate_flu_cpu,
        memory = abricate_flu_memory,
        docker = abricate_flu_docker,
        disk_size = abricate_flu_disk_size
    }
    # if flu_segment is not defined, assume the assembly is a full flu genome
    if (! defined(flu_segment)) {
      call extract_flu_segments_task.extract_flu_segments {
        input:
          assembly_fasta = assembly_fasta,
          flu_type = abricate_flu.abricate_flu_type,
          flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"])
      }
      call flu_antiviral.flu_antiviral_substitutions {
        input:
          na_segment_assembly = extract_flu_segments.seg_na_assembly,
          ha_segment_assembly = extract_flu_segments.seg_ha_assembly,
          pa_segment_assembly = extract_flu_segments.seg_pa_assembly,
          pb1_segment_assembly = extract_flu_segments.seg_pb1_assembly,
          pb2_segment_assembly = extract_flu_segments.seg_pb2_assembly,
          mp_segment_assembly = extract_flu_segments.seg_mp_assembly,
          abricate_flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]),
          irma_flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]),
          antiviral_aa_subs = antiviral_aa_subs,
          flu_h1_ha_ref = flu_h1_ha_ref,
          flu_h3_ha_ref = flu_h3_ha_ref,
          flu_n1_na_ref = flu_n1_na_ref,
          flu_n2_na_ref = flu_n2_na_ref,
          flu_pa_ref = flu_pa_ref,
          flu_pb1_ref = flu_pb1_ref,
          flu_pb2_ref = flu_pb2_ref,
          flu_h1n1_m2_ref = flu_h1n1_m2_ref,
          flu_h3n2_m2_ref = flu_h3n2_m2_ref
      }
    }
    # attempt to run nextclade for flu NA segment
    call set_organism_defaults.organism_parameters as set_flu_na_nextclade_values {
      input:
        organism = organism,
        flu_segment = "NA",
        flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"])
    }
    if (set_flu_na_nextclade_values.nextclade_dataset_tag != "NA") {
      call nextclade_task.nextclade_v3 as nextclade_flu_na {
        input:
          genome_fasta = select_first([extract_flu_segments.seg_na_assembly, assembly_fasta]),
          dataset_name = set_flu_na_nextclade_values.nextclade_dataset_name,
          dataset_tag = set_flu_na_nextclade_values.nextclade_dataset_tag,
          docker = nextclade_docker_image,
          cpu = nextclade_cpu,
          memory = nextclade_memory,
          disk_size = nextclade_disk_size
      }
      call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_na {
        input:
          nextclade_tsv = nextclade_flu_na.nextclade_tsv,
          organism = set_flu_na_nextclade_values.standardized_organism,
          docker = nextclade_output_parser_docker,
          cpu = nextclade_output_parser_cpu,
          memory = nextclade_output_parser_memory,
          disk_size = nextclade_output_parser_disk_size
      }
    }
    # attempt to run nextclade for flu HA segment
    call set_organism_defaults.organism_parameters as set_flu_ha_nextclade_values {
      input:
        organism = organism,
        flu_segment = "HA",
        flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"])
    }
    if (set_flu_ha_nextclade_values.nextclade_dataset_tag != "NA") {
      call nextclade_task.nextclade_v3 as nextclade_flu_ha {
        input:
          genome_fasta = select_first([extract_flu_segments.seg_ha_assembly, assembly_fasta]),
          dataset_name = set_flu_ha_nextclade_values.nextclade_dataset_name,
          dataset_tag = set_flu_ha_nextclade_values.nextclade_dataset_tag,
          docker = nextclade_docker_image,
          cpu = nextclade_cpu,
          memory = nextclade_memory,
          disk_size = nextclade_disk_size
      }
      call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_ha {
        input:
          nextclade_tsv = nextclade_flu_ha.nextclade_tsv,
          organism = set_flu_ha_nextclade_values.standardized_organism,
          docker = nextclade_output_parser_docker,
          cpu = nextclade_output_parser_cpu,
          memory = nextclade_output_parser_memory,
          disk_size = nextclade_output_parser_disk_size
      }
    }
    # only run GenoFLU and custom nextclade dataset if the subtype is H5N1 and the clade is 2.3.4.4b as they are specific to this subtype and clade.
    if (select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]) == "H5N1" && select_first([nextclade_output_parser_flu_ha.nextclade_clade, ""]) == "2.3.4.4b") {
      call genoflu_task.genoflu {
        input:
          assembly_fasta = assembly_fasta,
          samplename = samplename,
          min_percent_identity = genoflu_min_percent_identity,
          cross_reference = genoflu_cross_reference,
          cpu = genoflu_cpu,
          disk_size = genoflu_disk_size,
          docker = genoflu_docker,
          memory = genoflu_memory
      }
      call set_organism_defaults.organism_parameters as set_flu_h5n1_nextclade_values {
        input:
          organism = organism,
          flu_genoflu_genotype = genoflu.genoflu_genotype
      }
      if (genoflu.genoflu_genotype == "B3.13" || genoflu.genoflu_genotype == "D1.1" || defined(nextclade_custom_input_dataset)) {
        call nextclade_task.nextclade_v3 as nextclade_flu_h5n1 {
          input:
            genome_fasta = select_first([extract_flu_segments.concatenated_fasta, assembly_fasta]),
            custom_input_dataset = select_first([nextclade_custom_input_dataset, set_flu_h5n1_nextclade_values.nextclade_custom_dataset]),
            docker = nextclade_docker_image,
            cpu = nextclade_cpu,
            memory = nextclade_memory,
            disk_size = nextclade_disk_size
        }
        call nextclade_task.nextclade_output_parser as nextclade_output_parser_flu_h5n1 {
          input:
            nextclade_tsv = nextclade_flu_h5n1.nextclade_tsv,
            organism = set_flu_h5n1_nextclade_values.standardized_organism,
            docker = nextclade_output_parser_docker,
            cpu = nextclade_output_parser_cpu,
            memory = nextclade_output_parser_memory,
            disk_size = nextclade_output_parser_disk_size
        }
      }
    }
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      flu_segment = flu_segment,
      flu_subtype = select_first([flu_subtype, abricate_flu.abricate_flu_subtype, "N/A"]),
      reference_genome = reference_genome,
      genome_length_input = genome_length,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,
      vadr_max_length = vadr_max_length,
      vadr_skip_length = vadr_skip_length,
      vadr_options = vadr_opts,
      vadr_mem = vadr_memory
  }
  call consensus_qc_task.consensus_qc {
    input:
      assembly_fasta = assembly_fasta,
      reference_genome = organism_parameters.reference,
      genome_length = organism_parameters.genome_length
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta,
        docker = organism_parameters.pangolin_docker
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "measles") {
    if (organism_parameters.nextclade_dataset_tag != "NA") {
      call nextclade_task.nextclade_v3 {
        input:
          genome_fasta = select_first([extract_flu_segments.concatenated_fasta, assembly_fasta]),
          dataset_name = organism_parameters.nextclade_dataset_name,
          dataset_tag = organism_parameters.nextclade_dataset_tag,
          docker = nextclade_docker_image,
          cpu = nextclade_cpu,
          memory = nextclade_memory,
          disk_size = nextclade_disk_size
      }
      call nextclade_task.nextclade_output_parser {
        input:
          nextclade_tsv = nextclade_v3.nextclade_tsv,
          organism = organism_parameters.standardized_organism,
          docker = nextclade_output_parser_docker,
          cpu = nextclade_output_parser_cpu,
          memory = nextclade_output_parser_memory,
          disk_size = nextclade_output_parser_disk_size
      }
    }
  }
  # vadr task
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "flu") {
    call vadr_task.vadr {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = consensus_qc.number_ATCG,
        max_length = organism_parameters.vadr_maxlength,
        vadr_opts = organism_parameters.vadr_opts,
        skip_length = organism_parameters.vadr_skiplength,
        memory = organism_parameters.vadr_memory
    }
  }
  # QC check task
  if (defined(qc_check_table)) {
    call qc_check.qc_check_phb {
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
    # Extracted flu segments
    File? flu_segments_fasta_concatenated = extract_flu_segments.concatenated_fasta
    File? flu_ha_segment_fasta = extract_flu_segments.seg_ha_assembly
    File? flu_na_segment_fasta = extract_flu_segments.seg_na_assembly
    File? flu_pa_segment_fasta = extract_flu_segments.seg_pa_assembly
    File? flu_pb1_segment_fasta = extract_flu_segments.seg_pb1_assembly
    File? flu_pb2_segment_fasta = extract_flu_segments.seg_pb2_assembly
    File? flu_mp_segment_fasta = extract_flu_segments.seg_mp_assembly
    File? flu_np_segment_fasta = extract_flu_segments.seg_np_assembly
    File? flu_ns_segment_fasta = extract_flu_segments.seg_ns_assembly
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
    # Nextclade HA outputs
    File? nextclade_json_flu_ha = nextclade_flu_ha.nextclade_json
    File? auspice_json_flu_ha =  nextclade_flu_ha.auspice_json
    File? nextclade_tsv_flu_ha = nextclade_flu_ha.nextclade_tsv
    String? nextclade_ds_tag_flu_ha = set_flu_ha_nextclade_values.nextclade_dataset_tag
    String? nextclade_aa_subs_flu_ha = nextclade_output_parser_flu_ha.nextclade_aa_subs
    String? nextclade_aa_dels_flu_ha = nextclade_output_parser_flu_ha.nextclade_aa_dels
    String? nextclade_clade_flu_ha = nextclade_output_parser_flu_ha.nextclade_clade
    String? nextclade_qc_flu_ha = nextclade_output_parser_flu_ha.nextclade_qc
    # Nextclade NA outputs
    File? nextclade_json_flu_na = nextclade_flu_na.nextclade_json
    File? auspice_json_flu_na = nextclade_flu_na.auspice_json
    File? nextclade_tsv_flu_na = nextclade_flu_na.nextclade_tsv
    String? nextclade_ds_tag_flu_na = set_flu_na_nextclade_values.nextclade_dataset_tag
    String? nextclade_aa_subs_flu_na = nextclade_output_parser_flu_na.nextclade_aa_subs
    String? nextclade_aa_dels_flu_na = nextclade_output_parser_flu_na.nextclade_aa_dels
    String? nextclade_clade_flu_na = nextclade_output_parser_flu_na.nextclade_clade
    String? nextclade_qc_flu_na = nextclade_output_parser_flu_na.nextclade_qc
    # Nextclade H5N1 outputs
    File? nextclade_json_flu_h5n1 = nextclade_flu_h5n1.nextclade_json
    File? auspice_json_flu_h5n1 = nextclade_flu_h5n1.auspice_json
    File? nextclade_tsv_flu_h5n1 = nextclade_flu_h5n1.nextclade_tsv
    String? nextclade_aa_subs_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_aa_subs
    String? nextclade_aa_dels_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_aa_dels
    String? nextclade_clade_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_clade
    String? nextclade_qc_flu_h5n1 = nextclade_output_parser_flu_h5n1.nextclade_qc
    # VADR Annotation QC
    File?  vadr_alerts_list = vadr.alerts_list
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
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
    # GenoFLU outputs    
    String? genoflu_version = genoflu.genoflu_version
    String? genoflu_genotype = genoflu.genoflu_genotype
    String? genoflu_all_segments = genoflu.genoflu_all_segments
    File? genoflu_output_tsv = genoflu.genoflu_output_tsv
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = flu_antiviral_substitutions.flu_A_315675_resistance
    String? flu_amantadine_resistance = flu_antiviral_substitutions.flu_amantadine_resistance
    String? flu_compound_367_resistance = flu_antiviral_substitutions.flu_compound_367_resistance
    String? flu_favipiravir_resistance = flu_antiviral_substitutions.flu_favipiravir_resistance
    String? flu_fludase_resistance = flu_antiviral_substitutions.flu_fludase_resistance
    String? flu_L_742_001_resistance = flu_antiviral_substitutions.flu_L_742_001_resistance
    String? flu_laninamivir_resistance = flu_antiviral_substitutions.flu_laninamivir_resistance
    String? flu_peramivir_resistance = flu_antiviral_substitutions.flu_peramivir_resistance
    String? flu_pimodivir_resistance = flu_antiviral_substitutions.flu_pimodivir_resistance
    String? flu_rimantadine_resistance = flu_antiviral_substitutions.flu_rimantadine_resistance
    String? flu_oseltamivir_resistance = flu_antiviral_substitutions.flu_oseltamivir_resistance
    String? flu_xofluza_resistance = flu_antiviral_substitutions.flu_xofluza_resistance
    String? flu_zanamivir_resistance = flu_antiviral_substitutions.flu_zanamivir_resistance
  }
}