version 1.0

import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../utilities/wf_flu_track.wdl" as flu_track_wf
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../../tasks/species_typing/lentivirus/task_quasitools.wdl" as quasitools_task
import "../../tasks/quality_control/advanced_metrics/task_vadr.wdl" as vadr_task
import "../../tasks/quality_control/basic_statistics/task_gene_coverage.wdl" as gene_coverage_task

workflow morgana_magic {
  input {
    String samplename
    File? assembly_fasta
    String taxon_name
    String seq_method
    File read1
    File? read2
    Int? number_ATCG # needed for vadr
    # assembly metrics 
    Int? assembly_metrics_cpu
    Int? assembly_metrics_disk_size
    String? assembly_metrics_docker
    Int? assembly_metrics_memory
    # flu track - irma
    Int? irma_cpu
    Int? irma_disk_size
    String? irma_docker_image
    Boolean? irma_keep_ref_deletions
    Int? irma_memory
    # flu track - genoflu
    Int? genoflu_cpu
    File? genoflu_cross_reference
    Int? genoflu_disk_size
    String? genoflu_docker
    Int? genoflu_memory
    # flu track - abricate
    Int? abricate_flu_cpu
    Int? abricate_flu_disk_size
    String? abricate_flu_docker
    Int? abricate_flu_memory
    Int? abricate_flu_min_percent_coverage
    Int? abricate_flu_min_percent_identity
    Int? flu_track_min_depth
    # nextclade inputs
    String? nextclade_dataset_name
    String? nextclade_dataset_tag
    Int? nextclade_cpu
    Int? nextclade_disk_size
    String? nextclade_docker_image
    Int? nextclade_memory
    Int? nextclade_output_parser_cpu
    Int? nextclade_output_parser_disk_size
    String? nextclade_output_parser_docker
    Int? nextclade_output_parser_memory
    # pangolin inputs
    Int? pangolin_cpu
    Int? pangolin_disk_size
    String? pangolin_docker_image
    Int? pangolin_memory
    # gene coverage inputs
    File? reference_gene_locations_bed
    File? ivar_consensus_aligned_bam
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_skip_length
    String? vadr_options
    File? vadr_model_file
    Int? vadr_memory
    # Workflow Settings
    String? workflow_type
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = taxon_name,
      pangolin_docker_image = pangolin_docker_image
  }
  if (organism_parameters.standardized_organism == "flu") {
    if (workflow_type == "theiacov_fasta"){
      call vadr_task.vadr as vadr_fasta {
        input:
          genome_fasta = select_first([assembly_fasta]),
          assembly_length_unambiguous = select_first([number_ATCG]),
          vadr_opts = select_first([vadr_options, organism_parameters.vadr_opts]),
          vadr_model_file = select_first([vadr_model_file, organism_parameters.vadr_model_file]),
          max_length = select_first([vadr_max_length, organism_parameters.vadr_maxlength]),
          skip_length = select_first([vadr_skip_length, organism_parameters.vadr_skiplength]),
          memory = select_first([vadr_memory, organism_parameters.vadr_memory])
      }
      call flu_track_wf.flu_track as flu_track_fasta {
        input:
          samplename = samplename,
          assembly_fasta = select_first([assembly_fasta]),
          seq_method = seq_method,
          standardized_organism = organism_parameters.standardized_organism,
          assembly_metrics_cpu = assembly_metrics_cpu,
          assembly_metrics_disk_size = assembly_metrics_disk_size,
          assembly_metrics_docker = assembly_metrics_docker,
          assembly_metrics_memory = assembly_metrics_memory,
          irma_cpu = irma_cpu,
          irma_disk_size = irma_disk_size,
          irma_docker_image = irma_docker_image,        
          irma_keep_ref_deletions = irma_keep_ref_deletions,
          irma_memory = irma_memory,
          genoflu_cross_reference = genoflu_cross_reference,
          genoflu_cpu = genoflu_cpu,
          genoflu_disk_size = genoflu_disk_size,
          genoflu_docker = genoflu_docker,
          genoflu_memory = genoflu_memory,
          abricate_flu_cpu = abricate_flu_cpu,
          abricate_flu_disk_size = abricate_flu_disk_size,
          abricate_flu_docker = abricate_flu_docker,
          abricate_flu_memory = abricate_flu_memory,
          abricate_flu_min_percent_coverage = abricate_flu_min_percent_coverage,
          abricate_flu_min_percent_identity = abricate_flu_min_percent_identity,
          nextclade_cpu = nextclade_cpu,
          nextclade_disk_size = nextclade_disk_size,
          nextclade_docker_image = nextclade_docker_image,
          nextclade_memory = nextclade_memory,
          nextclade_output_parser_cpu = nextclade_output_parser_cpu,
          nextclade_output_parser_disk_size = nextclade_output_parser_disk_size,
          nextclade_output_parser_docker = nextclade_output_parser_docker,
          nextclade_output_parser_memory = nextclade_output_parser_memory,
          vadr_outputs_tgz = vadr.outputs_tgz
      }
    }
    if (workflow_type == "theiacov_pe") {
      call flu_track_wf.flu_track as flu_track_pe {
        input:
          samplename = samplename,
          read1 = read1,
          read2 = read2,
          seq_method = seq_method,
          standardized_organism = organism_parameters.standardized_organism,
          irma_min_consensus_support = select_first([flu_track_min_depth])
      }
      call vadr_task.vadr as vadr_pe {
        input:
          genome_fasta = select_first([flu_track_pe.irma_assembly_fasta_padded]),
          assembly_length_unambiguous = select_first([number_ATCG]),
          vadr_opts = select_first([vadr_options, organism_parameters.vadr_opts]),
          vadr_model_file = select_first([vadr_model_file, organism_parameters.vadr_model_file]),
          max_length = select_first([vadr_max_length, organism_parameters.vadr_maxlength]),
          skip_length = select_first([vadr_skip_length, organism_parameters.vadr_skiplength]),
          memory = select_first([vadr_memory, organism_parameters.vadr_memory])
      }
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles") {
    # tasks specific to MPXV, sars-cov-2, WNV, flu, rsv_a, and rsv_b
    call vadr_task.vadr  {
      input:
        genome_fasta = select_first([assembly_fasta]),
        assembly_length_unambiguous = select_first([number_ATCG]),
        vadr_opts = select_first([vadr_options, organism_parameters.vadr_opts]),
        vadr_model_file = select_first([vadr_model_file, organism_parameters.vadr_model_file]),
        max_length = select_first([vadr_max_length, organism_parameters.vadr_maxlength]),
        skip_length = select_first([vadr_skip_length, organism_parameters.vadr_skiplength]),
        memory = select_first([vadr_memory, organism_parameters.vadr_memory])
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || defined(reference_gene_locations_bed)) {
    # tasks specific to either sars-cov-2, MPXV, or any organism with a user-supplied reference gene locations bed file
    call gene_coverage_task.gene_coverage {
      input:
        bamfile = select_first([ivar_consensus_aligned_bam, flu_track_pe.irma_ha_bam, flu_track_pe.irma_na_bam, flu_track_fasta.irma_ha_bam, flu_track_fasta.irma_na_bam]),
        bedfile = select_first([reference_gene_locations_bed, organism_parameters.gene_locations_bed]),
        samplename = samplename,
        organism = select_first([organism_parameters.standardized_organism, taxon_name])
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = select_first([assembly_fasta]),
        docker = select_first([pangolin_docker_image, organism_parameters.pangolin_docker]),
        cpu = pangolin_cpu,
        disk_size = pangolin_disk_size,
        memory = pangolin_memory
    }
  }
  # Non-flu Nextclade
  if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "measles" ) {
    call nextclade_task.nextclade_v3 {
      input:
        genome_fasta = select_first([assembly_fasta]),
        dataset_name = select_first([nextclade_dataset_name, organism_parameters.nextclade_dataset_name]),
        dataset_tag = select_first([nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag]),
        cpu = nextclade_cpu,
        disk_size = nextclade_disk_size,
        docker = nextclade_docker_image,
        memory = nextclade_memory
    }
    call nextclade_task.nextclade_output_parser {
      input:
        nextclade_tsv = nextclade_v3.nextclade_tsv,
        organism = organism_parameters.standardized_organism,
        cpu = nextclade_output_parser_cpu,
        disk_size = nextclade_output_parser_disk_size,
        docker = nextclade_output_parser_docker,
        memory = nextclade_output_parser_memory
    }
  }
  if (organism_parameters.standardized_organism == "HIV") {
    call quasitools_task.quasitools {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename
    }
  }
  if (organism_parameters.standardized_organism == "rabies") {
    call nextclade_task.nextclade_v3_set as rabies_nextclade {
      input:
        genome_fastas = [select_first([assembly_fasta])],
        reference_tree_json = organism_parameters.nextclade_auspice_tree,
        gene_annotations_gff = organism_parameters.reference_gff,
        pathogen_json = organism_parameters.nextclade_pathogen_json,
        input_ref = organism_parameters.reference,
        cpu = nextclade_cpu,
        disk_size = nextclade_disk_size,
        docker = nextclade_docker_image,
        memory = nextclade_memory
    }
    call nextclade_task.nextclade_output_parser as rabies_output_parser {
      input:
        nextclade_tsv = rabies_nextclade.nextclade_tsv,
        organism = organism_parameters.standardized_organism,
        cpu = nextclade_output_parser_cpu,
        disk_size = nextclade_output_parser_disk_size,
        docker = nextclade_output_parser_docker,
        memory = nextclade_output_parser_memory
    }
  }
  output {
    String organism = select_first([organism_parameters.standardized_organism, taxon_name])
    # VADR outputs
    String? vadr_alerts_list = select_first([vadr.alerts_list, vadr_pe.alerts_list, vadr_fasta.alerts_list, ""])
    String? vadr_num_alerts = select_first([vadr.num_alerts, vadr_pe.num_alerts, vadr_fasta.num_alerts, ""])
    String? vadr_feature_tbl_pass = select_first([vadr.feature_tbl_pass, vadr_pe.feature_tbl_pass, vadr_fasta.feature_tbl_pass, ""])
    String? vadr_feature_tbl_fail = select_first([vadr.feature_tbl_fail, vadr_pe.feature_tbl_fail, vadr_fasta.feature_tbl_fail, ""])
    String? vadr_classification_summary_file = select_first([vadr.classification_summary_file, vadr_pe.classification_summary_file, vadr_fasta.classification_summary_file, ""])
    String? vadr_all_outputs_tar_gz = select_first([vadr.outputs_tgz, vadr_pe.outputs_tgz, vadr_fasta.outputs_tgz, ""])
    String? vadr_docker = select_first([vadr.vadr_docker, vadr_pe.vadr_docker, vadr_fasta.vadr_docker, ""])
    String? vadr_fastas_zip_archive = select_first([vadr.vadr_fastas_zip_archive, vadr_pe.vadr_fastas_zip_archive, vadr_fasta.vadr_fastas_zip_archive, ""])
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
    String nextclade_version = select_first([rabies_nextclade.nextclade_version, nextclade_v3.nextclade_version, flu_track_pe.nextclade_version, flu_track_fasta.nextclade_version, ""])
    String nextclade_docker = select_first([rabies_nextclade.nextclade_docker, nextclade_v3.nextclade_docker, flu_track_pe.nextclade_docker, flu_track_fasta.nextclade_docker, ""])
    String nextclade_ds_tag = select_first([nextclade_dataset_tag, organism_parameters.nextclade_dataset_tag, ""])
    # Nextclade outputs for non-flu
    File? nextclade_json = nextclade_v3.nextclade_json
    File? auspice_json = nextclade_v3.auspice_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    String? nextclade_aa_subs = nextclade_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels = nextclade_output_parser.nextclade_aa_dels
    String? nextclade_clade = nextclade_output_parser.nextclade_clade
    String? nextclade_lineage = nextclade_output_parser.nextclade_lineage
    String? nextclade_qc = nextclade_output_parser.nextclade_qc
    # Nextclade outputs for rabies
    File? nextclade_json_rabies = rabies_nextclade.nextclade_json
    File? auspice_json_rabies = rabies_nextclade.auspice_json
    File? nextclade_tsv_rabies = rabies_nextclade.nextclade_tsv
    String? nextclade_aa_subs_rabies = rabies_output_parser.nextclade_aa_subs
    String? nextclade_aa_dels_rabies = rabies_output_parser.nextclade_aa_dels
    String? nextclade_clade_rabies = rabies_output_parser.nextclade_clade
    String? nextclade_lineage_rabies = rabies_output_parser.nextclade_lineage
    String? nextclade_qc_rabies = rabies_output_parser.nextclade_qc
    # Nextclade outputs for flu H5N1
    String? nextclade_json_flu_h5n1 = select_first([flu_track_pe.nextclade_json_flu_h5n1, flu_track_fasta.nextclade_json_flu_h5n1, ""])
    String? auspice_json_flu_h5n1 = select_first([flu_track_pe.auspice_json_flu_h5n1, flu_track_fasta.auspice_json_flu_h5n1, ""])
    String? nextclade_tsv_flu_h5n1 = select_first([flu_track_pe.nextclade_tsv_flu_h5n1, flu_track_fasta.nextclade_tsv_flu_h5n1, ""])
    String? nextclade_aa_subs_flu_h5n1 = select_first([flu_track_pe.nextclade_aa_subs_flu_h5n1, flu_track_fasta.nextclade_aa_subs_flu_h5n1, ""])
    String? nextclade_aa_dels_flu_h5n1 = select_first([flu_track_pe.nextclade_aa_dels_flu_h5n1, flu_track_fasta.nextclade_aa_dels_flu_h5n1, ""])
    String? nextclade_clade_flu_h5n1 = select_first([flu_track_pe.nextclade_clade_flu_h5n1, flu_track_fasta.nextclade_clade_flu_h5n1, ""])
    String? nextclade_qc_flu_h5n1 = select_first([flu_track_pe.nextclade_qc_flu_h5n1, flu_track_fasta.nextclade_qc_flu_h5n1, ""])
    # Nextclade outputs for flu HA
    String? nextclade_json_flu_ha = select_first([flu_track_pe.nextclade_json_flu_ha, flu_track_fasta.nextclade_json_flu_ha, ""])
    String? auspice_json_flu_ha = select_first([flu_track_pe.auspice_json_flu_ha, flu_track_fasta.auspice_json_flu_ha, ""])
    String? nextclade_tsv_flu_ha = select_first([flu_track_pe.nextclade_tsv_flu_ha, flu_track_fasta.nextclade_tsv_flu_ha, ""])
    String? nextclade_ds_tag_flu_ha = select_first([flu_track_pe.nextclade_ds_tag_flu_ha, flu_track_fasta.nextclade_ds_tag_flu_ha, ""])
    String? nextclade_aa_subs_flu_ha = select_first([flu_track_pe.nextclade_aa_subs_flu_ha, flu_track_fasta.nextclade_aa_subs_flu_ha, ""])
    String? nextclade_aa_dels_flu_ha = select_first([flu_track_pe.nextclade_aa_dels_flu_ha, flu_track_fasta.nextclade_aa_dels_flu_ha, ""])
    String? nextclade_clade_flu_ha = select_first([flu_track_pe.nextclade_clade_flu_ha, flu_track_fasta.nextclade_clade_flu_ha, ""])
    String? nextclade_qc_flu_ha = select_first([flu_track_pe.nextclade_qc_flu_ha, flu_track_fasta.nextclade_qc_flu_ha, ""])
    # Nextclade outputs for flu NA
    String? nextclade_json_flu_na = select_first([flu_track_pe.nextclade_json_flu_na, flu_track_fasta.nextclade_json_flu_na, ""])
    String? auspice_json_flu_na = select_first([flu_track_pe.auspice_json_flu_na, flu_track_fasta.auspice_json_flu_na, ""])
    String? nextclade_tsv_flu_na = select_first([flu_track_pe.nextclade_tsv_flu_na, flu_track_fasta.nextclade_tsv_flu_na, ""])
    String? nextclade_ds_tag_flu_na = select_first([flu_track_pe.nextclade_ds_tag_flu_na, flu_track_fasta.nextclade_ds_tag_flu_na, ""])
    String? nextclade_aa_subs_flu_na = select_first([flu_track_pe.nextclade_aa_subs_flu_na, flu_track_fasta.nextclade_aa_subs_flu_na, ""])
    String? nextclade_aa_dels_flu_na = select_first([flu_track_pe.nextclade_aa_dels_flu_na, flu_track_fasta.nextclade_aa_dels_flu_na, ""])
    String? nextclade_clade_flu_na = select_first([flu_track_pe.nextclade_clade_flu_na, flu_track_fasta.nextclade_clade_flu_na, ""])
    String? nextclade_qc_flu_na = select_first([flu_track_pe.nextclade_qc_flu_na, flu_track_fasta.nextclade_qc_flu_na, ""])
    # Flu IRMA Outputs
    String? irma_version = select_first([flu_track_pe.irma_version, flu_track_fasta.irma_version, ""])
    String? irma_docker = select_first([flu_track_pe.irma_docker, flu_track_fasta.irma_docker, ""])
    String? irma_type = select_first([flu_track_pe.irma_type, flu_track_fasta.irma_type, ""])
    String? irma_subtype = select_first([flu_track_pe.irma_subtype, flu_track_fasta.irma_subtype, ""])
    String? irma_subtype_notes = select_first([flu_track_pe.irma_subtype_notes, flu_track_fasta.irma_subtype_notes, ""])
    String? irma_assembly_fasta = select_first([flu_track_pe.irma_assembly_fasta, flu_track_fasta.irma_assembly_fasta, ""])
    String? flu_assembly_fasta_concatenated = select_first([flu_track_pe.flu_assembly_fasta_concatenated, flu_track_fasta.flu_assembly_fasta_concatenated, ""])
    String? irma_minimum_consensus_support = select_first([flu_track_pe.irma_minimum_consensus_support, flu_track_fasta.irma_minimum_consensus_support, ""])
    String? ha_na_assembly_coverage = select_first([flu_track_pe.ha_na_assembly_coverage, flu_track_fasta.ha_na_assembly_coverage, ""])
    String? flu_ha_segment_fasta = select_first([flu_track_pe.flu_ha_segment_fasta, flu_track_fasta.flu_ha_segment_fasta, ""])
    String? flu_na_segment_fasta = select_first([flu_track_pe.flu_na_segment_fasta, flu_track_fasta.flu_na_segment_fasta, ""])
    String? flu_pa_segment_fasta = select_first([flu_track_pe.flu_pa_segment_fasta, flu_track_fasta.flu_pa_segment_fasta, ""])
    String? flu_pb1_segment_fasta = select_first([flu_track_pe.flu_pb1_segment_fasta, flu_track_fasta.flu_pb1_segment_fasta, ""])
    String? flu_pb2_segment_fasta = select_first([flu_track_pe.flu_pb2_segment_fasta, flu_track_fasta.flu_pb2_segment_fasta, ""])
    String? flu_mp_segment_fasta = select_first([flu_track_pe.flu_mp_segment_fasta, flu_track_fasta.flu_mp_segment_fasta, ""])
    String? flu_np_segment_fasta = select_first([flu_track_pe.flu_np_segment_fasta, flu_track_fasta.flu_np_segment_fasta, ""])
    String? flu_ns_segment_fasta = select_first([flu_track_pe.flu_ns_segment_fasta, flu_track_fasta.flu_ns_segment_fasta, ""])
    # Flu GenoFLU Outputs
    String? genoflu_version = select_first([flu_track_pe.genoflu_version, flu_track_fasta.genoflu_version, ""])
    String? genoflu_genotype = select_first([flu_track_pe.genoflu_genotype, flu_track_fasta.genoflu_genotype, ""])
    String? genoflu_all_segments = select_first([flu_track_pe.genoflu_all_segments, flu_track_fasta.genoflu_all_segments, ""])
    String? genoflu_output_tsv = select_first([flu_track_pe.genoflu_output_tsv, flu_track_fasta.genoflu_output_tsv, ""])
    # Flu Abricate Outputs
    String? abricate_flu_type = select_first([flu_track_pe.abricate_flu_type, flu_track_fasta.abricate_flu_type, ""])
    String? abricate_flu_subtype =  select_first([flu_track_pe.abricate_flu_subtype, flu_track_fasta.abricate_flu_subtype, ""])
    String? abricate_flu_results = select_first([flu_track_pe.abricate_flu_results, flu_track_fasta.abricate_flu_results, ""])
    String? abricate_flu_database =  select_first([flu_track_pe.abricate_flu_database, flu_track_fasta.abricate_flu_database, ""])
    String? abricate_flu_version = select_first([flu_track_pe.abricate_flu_version, flu_track_fasta.abricate_flu_version, ""])
    # Flu Antiviral Substitution Outputs
    String? flu_A_315675_resistance = select_first([flu_track_pe.flu_A_315675_resistance, flu_track_fasta.flu_A_315675_resistance, ""])
    String? flu_amantadine_resistance = select_first([flu_track_pe.flu_amantadine_resistance, flu_track_fasta.flu_amantadine_resistance, ""])
    String? flu_compound_367_resistance = select_first([flu_track_pe.flu_compound_367_resistance, flu_track_fasta.flu_compound_367_resistance, ""])
    String? flu_favipiravir_resistance = select_first([flu_track_pe.flu_favipiravir_resistance, flu_track_fasta.flu_favipiravir_resistance, ""])
    String? flu_fludase_resistance = select_first([flu_track_pe.flu_fludase_resistance, flu_track_fasta.flu_fludase_resistance, ""])
    String? flu_L_742_001_resistance = select_first([flu_track_pe.flu_L_742_001_resistance, flu_track_fasta.flu_L_742_001_resistance, ""])
    String? flu_laninamivir_resistance = select_first([flu_track_pe.flu_laninamivir_resistance, flu_track_fasta.flu_laninamivir_resistance, ""])
    String? flu_peramivir_resistance = select_first([flu_track_pe.flu_peramivir_resistance, flu_track_fasta.flu_peramivir_resistance, ""])
    String? flu_pimodivir_resistance = select_first([flu_track_pe.flu_pimodivir_resistance, flu_track_fasta.flu_pimodivir_resistance, ""])
    String? flu_rimantadine_resistance = select_first([flu_track_pe.flu_rimantadine_resistance, flu_track_fasta.flu_rimantadine_resistance, ""])
    String? flu_oseltamivir_resistance = select_first([flu_track_pe.flu_oseltamivir_resistance, flu_track_fasta.flu_oseltamivir_resistance, ""])
    String? flu_xofluza_resistance = select_first([flu_track_pe.flu_xofluza_resistance, flu_track_fasta.flu_xofluza_resistance, ""])
    String? flu_zanamivir_resistance = select_first([flu_track_pe.flu_zanamivir_resistance, flu_track_fasta.flu_zanamivir_resistance, ""])
    # HIV Quasitools Outputs
    String? quasitools_version = quasitools.quasitools_version
    String? quasitools_date = quasitools.quasitools_date
    File? quasitools_coverage_file = quasitools.coverage_file
    File? quasitools_dr_report = quasitools.dr_report
    File? quasitools_hydra_vcf = quasitools.hydra_vcf
    File? quasitools_mutations_report = quasitools.mutations_report
    # Gene Coverage Outputs
    Float? sc2_s_gene_mean_coverage = gene_coverage.sc2_s_gene_depth
    Float? sc2_s_gene_percent_coverage = gene_coverage.sc2_s_gene_percent_coverage
    File? est_percent_gene_coverage_tsv = gene_coverage.est_percent_gene_coverage_tsv
    # Flu Track Percentage Mapped Reads
    String? percentage_mapped_reads = select_first([flu_track_pe.percentage_mapped_reads, flu_track_fasta.percentage_mapped_reads, ""])
  }
}