version 1.0

task czgenepi_wrangling {
  input {
    File full_terra_table
    Array[String] sample_names
    String? terra_project
    String? terra_workspace
    String? terra_table
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
    echo "Sample Name (from FASTA),Private ID,GISAID ID (Public ID) - Optional,GenBank Accession (Public ID) - Optional,Collection Date,Collection Location,Sequencing Date - Optional,Sample is Private"

    # parse terra table for data
    python3 <<CODE
    import pandas as pd
    import numpy as np
    import re

    # read in terra table
    table = pd.read_csv("~{full_terra_table}", delimiter='\t', header=0, dtype={"~{sample_names}": 'str'}) # ensure sample_id is always a string



    CODE
  >>>
  output {
    Array[String] fastas = read_lines("file_list.txt")
    File concatenated_metadata = "file"
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