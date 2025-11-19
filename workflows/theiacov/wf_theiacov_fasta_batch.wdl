version 1.0

import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_handling/task_theiacov_fasta_batch.wdl" as theiacov_fasta_wrangling_task
import "../../tasks/utilities/file_handling/task_cat_files.wdl" as concatenate
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults
import "../utilities/wf_morgana_magic.wdl" as morgana_magic_wf

workflow theiacov_fasta_batch {
  meta {
    description: "TheiaCoV_FASTA for multiple samples"
  }
  input {
    Array[String] samplenames
    Array[File] assembly_fastas
    String organism = "sars-cov-2"
    # nextclade inputs
    String? nextclade_dataset_tag
    String? nextclade_dataset_name
    # pangolin inputs
    String? pangolin_docker
    # workspace values
    String table_name
    String workspace_name
    String project_name
    String bucket_name
  }
  call set_organism_defaults.organism_parameters {
    input:
      organism = organism,
      nextclade_dataset_tag_input = nextclade_dataset_tag,
      nextclade_dataset_name_input = nextclade_dataset_name,
      pangolin_docker_image = pangolin_docker
  }
  call concatenate.cat_files_fasta {
    input: 
      files_to_cat = assembly_fastas,
      headers = samplenames,
      concatenated_file_name = "concatenated_assemblies.fasta"
  }
  call morgana_magic_wf.morgana_magic {
    input:
      samplename = "concatenated_assemblies",
      assembly_fasta = cat_files_fasta.concatenated_files,
      taxon_name = organism_parameters.standardized_organism,
      nextclade_dataset_name = organism_parameters.nextclade_dataset_name,
      nextclade_dataset_tag = organism_parameters.nextclade_dataset_tag,
      pangolin_docker_image = organism_parameters.pangolin_docker,
      seq_method = "NA",
      # Setting flu_track related inputs to default values as they are not utilized in TheiaCov, decreasing external input bloat
      assembly_metrics_cpu = 0,
      assembly_metrics_disk_size = 0,
      assembly_metrics_docker = "",
      assembly_metrics_memory = 0,
      irma_cpu = 0,
      irma_disk_size = 0,
      irma_docker_image = "",        
      irma_keep_ref_deletions = false,
      irma_memory = 0,
      genoflu_cross_reference = "",
      genoflu_cpu = 0,
      genoflu_disk_size = 0,
      genoflu_docker = "",
      genoflu_memory = 0,
      abricate_flu_cpu = 0,
      abricate_flu_disk_size = 0,
      abricate_flu_docker = "",
      abricate_flu_memory = 0,
      abricate_flu_min_percent_coverage = 0,
      abricate_flu_min_percent_identity = 0,
      flu_track_antiviral_aa_subs = "",
      nextclade_cpu = 0,
      nextclade_disk_size = 0,
      nextclade_docker_image = "",
      nextclade_memory = 0,
      nextclade_custom_input_dataset = "",
      nextclade_output_parser_cpu = 0,
      nextclade_output_parser_disk_size = 0,
      nextclade_output_parser_docker = "",
      nextclade_output_parser_memory = 0,
      workflow_type = "theiacov_fasta_batch"
  }
  call versioning.version_capture {
    input:
  }
  call theiacov_fasta_wrangling_task.sm_theiacov_fasta_wrangling {
    input:
      table_name = table_name,
      workspace_name = workspace_name,
      project_name = project_name,
      bucket_name = bucket_name,
      samplenames = samplenames,
      organism = organism_parameters.standardized_organism,
      nextclade_tsv = morgana_magic.nextclade_tsv,
      nextclade_docker = morgana_magic.nextclade_docker,
      nextclade_version = morgana_magic.nextclade_version,
      nextclade_ds_tag = organism_parameters.nextclade_dataset_tag,
      nextclade_json = morgana_magic.nextclade_json,
      pango_lineage_report = morgana_magic.pango_lineage_report,
      pangolin_docker = morgana_magic.pangolin_docker,
      theiacov_fasta_analysis_date = version_capture.date,
      theiacov_fasta_version = version_capture.phb_version
  }
  output {
    # Version Capture
    String theiacov_fasta_batch_version = version_capture.phb_version
    String theiacov_fasta_batch_analysis_date = version_capture.date
    # Pangolin outputs
    File? pango_lineage_report = morgana_magic.pango_lineage_report
    # Nextclade outputs
    File? nextclade_json = morgana_magic.nextclade_json
    File? nextclade_tsv = morgana_magic.nextclade_tsv
    # Wrangling outputs
    File datatable = sm_theiacov_fasta_wrangling.terra_table
  }
}