version 1.0

task czgenepi_wrangling {
  input {
    File full_terra_table
    Array[String] sample_names
    String assembly_fasta_column_name
    String collection_date_column_name
    String private_id_column_name

    # collection location
    String continent_column_name
    String country_column_name
    String state_column_name
    String county_column_name

    # optional inputs
    String gisaid_id_column_name
    String genbank_accession_column_name
    String sequencing_date_column_name
    String sample_is_private_column_name

    # runtime
    Int disk_size = 100
  }
  command <<<
    # create output file header
    echo "Sample Name (from FASTA),Private ID,GISAID ID (Public ID) - Optional,GenBank Accession (Public ID) - Optional,Collection Date,Collection Location,Sequencing Date - Optional,Sample is Private" > czgenepi_prep_metadata.csv

    # parse terra table for data
    python3 <<CODE
    import pandas as pd
    import numpy as np
    import re

    # read in terra table and set the private_id_column_name (defaulted to the terra table id) column to be only strings
    table = pd.read_csv("~{full_terra_table}", delimiter='\t', header=0, dtype={"~{private_id_column_name}": 'str'})

    # extract the samples for upload from the entire table
    table = table[table["~{private_id_column_name}"].isin("~{sep='*' sample_names}".split("*"))]

    # set all column headers to lowercase
    table.columns = table.columns.str.lower()

    # make lists for the columns needed for the metadata spreadsheet
    REQUIRED_COLUMNS = ["~{collection_date_column_name}", "~{private_id_column_name}", "~{continent_column_name}", "~{country_column_name}", "~{state_column_name}"]
    OPTIONAL_COLUMNS = ["~{county_column_name}", "~{gisaid_id_column_name}", "~{genbank_accession_column_name}", "~{sequencing_date_column_name}", "~{sample_is_private_column_name}"]

    metadata = table[REQUIRED_COLUMNS].copy()
    for column in OPTIONAL_COLUMNS:
      if column in table.columns:
        metadata[column] = table[column]
      else:
        if column == "sample_is_private":
          # by default, all samples will be set as private
          metadata[column] = "Yes" 
        else:
          metadata[column] = ""

    # write the output to a csv file
    metadata.to_csv("czgenepi_prep_metadata.csv", mode="a", index=False)
    
    # create a list of the assembly fastas for concatenation
    table.to_csv("file_list.txt", columns=["~{assembly_fasta_column_name}"], index=False, header=False)

    CODE
  >>>
  output {
    Array[String] fastas = read_lines("file_list.txt")
    File concatenated_metadata = "czgenepi_prep_metadata.csv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-08-08-2"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}