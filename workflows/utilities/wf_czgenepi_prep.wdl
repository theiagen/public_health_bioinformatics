version 1.0

import "../../tasks/utilities/task_czgenepi_wrangling.wdl" as czgenepi_wrangling_task
import "../../tasks/utilities/task_download_terra_table.wdl" as download_table
import "../../tasks/task_versioning.wdl" as versioning

workflow czgenepi_prep {
  input {
    Array[String] sample_names

    # downloading table information
    String terra_project_name
    String terra_workspace_name
    String terra_table_name
    
    # required columns
    String assembly_fasta_column_name = "assembly_fasta"
    String collection_date_column_name = "collection_date"
    String private_id_column_name = terra_table_name + "_id"

    # collection location - required
    String continent_column_name = "continent"
    String country_column_name = "country"
    String state_column_name = "state"
    String county_column_name = "county"

    # optional columns
    String gisaid_id_column_name = "gisaid_accession"
    String genbank_accession_column_name = "genbank_accession"
    String sequencing_date_column_name = "sequencing_date"
    String sample_is_private_column_name = "sample_is_private"
  }
  call versioning.version_capture{
    input:
  }
  call download_table.download_terra_table {
    input:
      terra_project_name = terra_project_name,
      terra_workspace_name = terra_workspace_name,
      terra_table_name = terra_table_name
  }
  call czgenepi_wrangling_task.czgenepi_wrangling {
    input:
      full_terra_table = download_terra_table.terra_table,
      sample_names = sample_names,
      terra_table_name = terra_table_name,
      assembly_fasta_column_name = assembly_fasta_column_name,
      collection_date_column_name = collection_date_column_name,
      private_id_column_name = private_id_column_name,
      continent_column_name = continent_column_name,
      country_column_name = country_column_name,
      state_column_name = state_column_name,
      county_column_name = county_column_name,
      gisaid_id_column_name = gisaid_id_column_name,
      genbank_accession_column_name = genbank_accession_column_name,
      sequencing_date_column_name = sequencing_date_column_name,
      sample_is_private_column_name = sample_is_private_column_name
  }
  output {
    File concatenated_czgenepi_fasta = czgenepi_wrangling.concatenated_fasta
    File concatenated_czgenepi_metadata = czgenepi_wrangling.concatenated_metadata
    String czgenepi_prep_version = version_capture.phb_version
    String czgenepi_prep_analysis_date = version_capture.date
  }
}