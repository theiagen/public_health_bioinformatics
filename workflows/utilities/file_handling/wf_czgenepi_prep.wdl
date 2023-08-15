version 1.0

import "../../../tasks/utilities/file_handling/task_czgenepi_wrangling.wdl" as czgenepi_wrangling_task


workflow czgenepi_prep {
  input {
    Array[String] sample_names
    Array[File] assembly_fasta
    Array[String] collection_date
    Array[String] private_id

    # collection location
    Array[String] continent
    Array[String] country
    Array[String] state
    Array[String]? county

    # optional inputs
    Array[String]? gisaid_virus_name
    Array[String]? sequencing_date
    Array[String]? sample_is_private
  }
  call czgenepi_wrangling_task.czgenepi_wrangling {
    input:
      sample_names = sample_names,
      assembly_fasta = assembly_fasta,
      collection_date = collection_date,
      private_id = private_id,
      continent = continent,
      country = country,
      state = state,
      county = county,
      gisaid_virus_name = gisaid_virus_name,
      sequencing_date = sequencing_date,
      sample_is_private = sample_is_private
  }
  output {
    File concatenated_fasta = czgenepi_wrangling.concatenated_fasta
    File concatenated_metadata = czgenepi_wrangling.concatenated_metadata
  }
}