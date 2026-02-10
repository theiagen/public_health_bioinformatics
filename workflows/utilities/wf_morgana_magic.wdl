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
    File assembly_fasta
    String taxon_name
    String seq_method
    File? read1
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
    # flu track - antiviral substitutions
    String? flu_track_antiviral_aa_subs
    # nextclade inputs
    String? nextclade_dataset_name
    String? nextclade_dataset_tag
    File? nextclade_custom_input_dataset
    File? nextclade_reference_gff_file
    File? nextclade_auspice_reference_tree_json
    File? nextclade_pathogen_json
    File? nextclade_input_ref
    String? nextclade_verbosity
    Int? nextclade_cpu
    Int? nextclade_disk_size
    String? nextclade_docker_image
    Int? nextclade_memory
    Int? nextclade_output_parser_cpu
    Int? nextclade_output_parser_disk_size
    String? nextclade_output_parser_docker
    Int? nextclade_output_parser_memory
    # pangolin inputs
    String? pangolin_analysis_mode
    Boolean? pangolin_expanded_lineage
    Float? pangolin_max_ambig
    Int? pangolin_min_length
    String? pangolin_arguments
    Boolean? pangolin_skip_designation_cache
    Boolean? pangolin_skip_scorpio
    Int? pangolin_cpu
    Int? pangolin_disk_size
    String? pangolin_docker_image
    Int? pangolin_memory
    # gene coverage inputs
    File? reference_gene_locations_bed
    File? gene_coverage_bam
    Int? gene_coverage_min_depth
    Int? sc2_s_gene_start
    Int? sc2_s_gene_stop
    Int? gene_coverage_cpu
    Int? gene_coverage_disk_size
    Int? gene_coverage_memory
    String? gene_coverage_docker
    # quasitools inputs
    Int? quasitools_cpu
    Int? quasitools_memory
    Int? quasitools_disk_size
    String? quasitools_docker
    # vadr parameters
    Int? vadr_max_length
    Int? vadr_min_length
    Int? vadr_skip_length
    String? vadr_options
    File? vadr_model_file
    Int? vadr_memory
    Int? vadr_cpu
    Int? vadr_disk_size
    String? vadr_docker_image
    # Workflow Settings
    String workflow_type
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = taxon_name,
      pangolin_docker_image = pangolin_docker_image,
      gene_locations_bed_file = reference_gene_locations_bed,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,     
      vadr_max_length = vadr_max_length,
      vadr_skip_length = vadr_skip_length,
      vadr_options = vadr_options,
      vadr_model = vadr_model_file,
      vadr_mem = vadr_memory,
  }
  if (workflow_type != "theiacov_fasta_batch" && (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "flu" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "WNV" || organism_parameters.standardized_organism == "mumps" || organism_parameters.standardized_organism == "rubella" || organism_parameters.standardized_organism == "measles")) {
    # tasks specific to MPXV, sars-cov-2, WNV, flu, rsv_a, and rsv_b
    call vadr_task.vadr  {
      input:
        genome_fasta = assembly_fasta,
        assembly_length_unambiguous = select_first([number_ATCG]),
        vadr_opts = organism_parameters.vadr_opts,
        vadr_model_file = organism_parameters.vadr_model_file,
        max_length = organism_parameters.vadr_maxlength,
        min_length = vadr_min_length,
        skip_length = organism_parameters.vadr_skiplength,
        memory = organism_parameters.vadr_memory,
        cpu = vadr_cpu,
        disk_size = vadr_disk_size,
        docker = vadr_docker_image
    }
  }
  if (organism_parameters.standardized_organism == "flu" && (workflow_type == "theiaviral" || workflow_type == "theiaviral_panel")) {
    call flu_track_wf.flu_track {
      input:
        samplename = samplename,
        assembly_fasta = assembly_fasta,
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
        antiviral_aa_subs = flu_track_antiviral_aa_subs,
        nextclade_cpu = nextclade_cpu,
        nextclade_disk_size = nextclade_disk_size,
        nextclade_docker_image = nextclade_docker_image,
        nextclade_memory = nextclade_memory,
        nextclade_custom_input_dataset = nextclade_custom_input_dataset,
        nextclade_output_parser_cpu = nextclade_output_parser_cpu,
        nextclade_output_parser_disk_size = nextclade_output_parser_disk_size,
        nextclade_output_parser_docker = nextclade_output_parser_docker,
        nextclade_output_parser_memory = nextclade_output_parser_memory,
        vadr_outputs_tgz = vadr.outputs_tgz
    }
  }
  if ((workflow_type == "theiacov_pe" || workflow_type == "theiacov_se" || workflow_type == "theiacov_ont" || workflow_type == "theiacov_clearlabs") && (organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "MPXV" || defined(organism_parameters.gene_locations_bed))) {
    if (basename(organism_parameters.gene_locations_bed) != "empty.bed") {
      # tasks specific to either sars-cov-2, MPXV, or any organism with a user-supplied reference gene locations bed file
      call gene_coverage_task.gene_coverage {
        input:
          bamfile = select_first([gene_coverage_bam]),
          bedfile = organism_parameters.gene_locations_bed,
          samplename = samplename,
          organism = organism_parameters.standardized_organism,
          sc2_s_gene_start = sc2_s_gene_start,
          sc2_s_gene_stop = sc2_s_gene_stop,
          min_depth = gene_coverage_min_depth,
          cpu = gene_coverage_cpu,
          disk_size = gene_coverage_disk_size,
          docker = gene_coverage_docker,
          memory = gene_coverage_memory
      }
    }
  }
  if (organism_parameters.standardized_organism == "sars-cov-2") {
    call pangolin.pangolin4 {
      input:
        samplename = samplename,
        fasta = assembly_fasta,
        analysis_mode = pangolin_analysis_mode,
        expanded_lineage = pangolin_expanded_lineage,
        max_ambig = pangolin_max_ambig,
        min_length = pangolin_min_length,
        skip_designation_cache = pangolin_skip_designation_cache,
        skip_scorpio = pangolin_skip_scorpio,
        pangolin_arguments = pangolin_arguments,
        docker = organism_parameters.pangolin_docker,
        cpu = pangolin_cpu,
        disk_size = pangolin_disk_size,
        memory = pangolin_memory
    }
  }
  # Non-flu Nextclade
  if (organism_parameters.standardized_organism == "MPXV" || organism_parameters.standardized_organism == "sars-cov-2" || organism_parameters.standardized_organism == "rsv_a" || organism_parameters.standardized_organism == "rsv_b" || organism_parameters.standardized_organism == "measles" ) {
    call nextclade_task.nextclade_v3 {
      input:
        genome_fasta = assembly_fasta,
        dataset_name = organism_parameters.nextclade_dataset_name,
        dataset_tag = organism_parameters.nextclade_dataset_tag,
        auspice_reference_tree_json = nextclade_auspice_reference_tree_json,
        gene_annotations_gff = nextclade_reference_gff_file,
        nextclade_pathogen_json = nextclade_pathogen_json, 
        input_ref = nextclade_input_ref,
        verbosity = nextclade_verbosity,
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
  if (defined(read1) && organism_parameters.standardized_organism == "HIV") {
    call quasitools_task.quasitools {
      input:
        read1 = select_first([read1]),
        read2 = read2,
        samplename = samplename,
        cpu = quasitools_cpu,
        memory = quasitools_memory,
        disk_size = quasitools_disk_size,
        docker = quasitools_docker
    }
  }
  if (organism_parameters.standardized_organism == "rabies") {
    call nextclade_task.nextclade_v3_set as rabies_nextclade {
      input:
        genome_fastas = [assembly_fasta],
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
    String organism = organism_parameters.standardized_organism
    # VADR outputs
    File? vadr_alerts_list = vadr.alerts_list
    String? vadr_num_alerts = vadr.num_alerts
    File? vadr_feature_tbl_pass = vadr.feature_tbl_pass
    File? vadr_feature_tbl_fail = vadr.feature_tbl_fail
    File? vadr_classification_summary_file = vadr.classification_summary_file
    File? vadr_all_outputs_tar_gz = vadr.outputs_tgz
    String? vadr_docker = vadr.vadr_docker
    File? vadr_fastas_zip_archive = vadr.vadr_fastas_zip_archive
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
    String nextclade_version = select_first([rabies_nextclade.nextclade_version, nextclade_v3.nextclade_version, flu_track.nextclade_version, ""])
    String nextclade_docker = select_first([rabies_nextclade.nextclade_docker, nextclade_v3.nextclade_docker, flu_track.nextclade_docker, ""])
    String nextclade_ds_tag = organism_parameters.nextclade_dataset_tag
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
    # Flu IRMA Outputs
    String? irma_version = flu_track.irma_version
    String? irma_docker = flu_track.irma_docker
    String? irma_type = flu_track.irma_type
    String? irma_subtype = flu_track.irma_subtype
    String? irma_subtype_notes = flu_track.irma_subtype_notes
    File? irma_assembly_fasta = flu_track.irma_assembly_fasta
    File? irma_aligned_fastqs = flu_track.irma_aligned_fastqs
    File? irma_read1_aligned = flu_track.irma_read1_aligned
    File? irma_read2_aligned = flu_track.irma_read2_aligned
    File? flu_assembly_fasta_concatenated = flu_track.flu_assembly_fasta_concatenated
    Int? irma_minimum_consensus_support = flu_track.irma_minimum_consensus_support
    String? ha_na_assembly_coverage = flu_track.ha_na_assembly_coverage
    File? flu_ha_segment_fasta = flu_track.flu_ha_segment_fasta
    File? flu_na_segment_fasta = flu_track.flu_na_segment_fasta
    File? flu_pa_segment_fasta = flu_track.flu_pa_segment_fasta
    File? flu_pb1_segment_fasta = flu_track.flu_pb1_segment_fasta
    File? flu_pb2_segment_fasta = flu_track.flu_pb2_segment_fasta
    File? flu_mp_segment_fasta = flu_track.flu_mp_segment_fasta
    File? flu_np_segment_fasta = flu_track.flu_np_segment_fasta
    File? flu_ns_segment_fasta = flu_track.flu_ns_segment_fasta
    # Flu GenoFLU Outputs
    String? genoflu_version = flu_track.genoflu_version
    String? genoflu_genotype = flu_track.genoflu_genotype
    String? genoflu_all_segments = flu_track.genoflu_all_segments
    File? genoflu_output_tsv = flu_track.genoflu_output_tsv
    # Flu Abricate Outputs
    String? abricate_flu_type = flu_track.abricate_flu_type
    String? abricate_flu_subtype = flu_track.abricate_flu_subtype
    File? abricate_flu_results = flu_track.abricate_flu_results
    String? abricate_flu_database = flu_track.abricate_flu_database
    String? abricate_flu_version = flu_track.abricate_flu_version
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
    String? percentage_mapped_reads = flu_track.percentage_mapped_reads
  }
}