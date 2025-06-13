version 1.0

import "../../tasks/species_typing/betacoronavirus/task_pangolin.wdl" as pangolin_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/taxon_id/task_nextclade.wdl" as nextclade_task
import "../../tasks/utilities/data_handling/task_theiacov_fasta_batch.wdl" as theiacov_fasta_wrangling_task
import "../../tasks/utilities/file_handling/task_cat_files.wdl" as concatenate
import "../utilities/wf_organism_parameters.wdl" as set_organism_defaults

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
  if (organism == "sars-cov-2") {
    # sars-cov-2 specific tasks
    call pangolin_task.pangolin4 {
      input:
        samplename = "concatenated_assemblies",
        fasta = cat_files_fasta.concatenated_files,
        docker = organism_parameters.pangolin_docker
    }
  }
  if (organism == "MPXV" || organism == "sars-cov-2" || organism_parameters.standardized_organism == "measles") {
    # tasks specific to either MPXV or sars-cov-2 
    call nextclade_task.nextclade_v3 {
      input:
      genome_fasta = cat_files_fasta.concatenated_files,
      dataset_name = organism_parameters.nextclade_dataset_name,
      dataset_tag = organism_parameters.nextclade_dataset_tag
    }
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
      nextclade_tsv = nextclade_v3.nextclade_tsv,
      nextclade_docker = nextclade_v3.nextclade_docker,
      nextclade_version = nextclade_v3.nextclade_version,
      nextclade_ds_tag = organism_parameters.nextclade_dataset_tag,
      nextclade_json = nextclade_v3.nextclade_json,
      pango_lineage_report = pangolin4.pango_lineage_report,
      pangolin_docker = pangolin4.pangolin_docker,
      theiacov_fasta_analysis_date = version_capture.date,
      theiacov_fasta_version = version_capture.phb_version
  }
  output {
    # Version Capture
    String theiacov_fasta_batch_version = version_capture.phb_version
    String theiacov_fasta_batch_analysis_date = version_capture.date
    # Pangolin outputs
    File? pango_lineage_report = pangolin4.pango_lineage_report
    # Nextclade outputs
    File? nextclade_json = nextclade_v3.nextclade_json
    File? nextclade_tsv = nextclade_v3.nextclade_tsv
    # Wrangling outputs
    File datatable = sm_theiacov_fasta_wrangling.terra_table
  }
}