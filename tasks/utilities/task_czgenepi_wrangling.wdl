version 1.0

task czgenepi_wrangling {
  input {
    File full_terra_table
    String terra_table_name

    # required fields
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
    String genbank_accession_column_name
    String sequencing_date_column_name

    # setting columns
    String organism
    Boolean is_private = true

    # runtime
    Int disk_size = 100
  }
  command <<<
    export is_private=~{is_private}
   
    # parse terra table for data
    python3 <<CODE
    import pandas as pd
    import numpy as np
    import re
    import os

    # set a function to grab only the year from the date
    def year_getter(date):
      r = re.compile('^\d{4}-\d{2}-\d{2}')      
      if pd.isna(date):
        print("Incorrect collection date format; collection date must be in YYYY-MM-DD format. Invalid date was: NaN")
      elif r.match(date) is None:
        print("Incorrect collection date format; collection date must be in YYYY-MM-DD format. Invalid date was: " + str(date))
        return np.nan
      else:
        return date.split("-")[0]

    if os.environ["is_private"] == "true":
      is_private = "Yes"
    else:
      is_private = "No"

    # make lists for the columns needed for the metadata spreadsheet
    REQUIRED_COLUMNS = ["~{terra_table_name}_id", "~{collection_date_column_name}", "~{private_id_column_name}", "~{continent_column_name}", "~{country_column_name}", "~{state_column_name}"]
    OPTIONAL_COLUMNS_SC2 = ["~{county_column_name}", "~{sequencing_date_column_name}"]
    OPTIONAL_COLUMNS_MPOX = ["~{county_column_name}", "~{genbank_accession_column_name}", "~{sequencing_date_column_name}"]
    OUTPUT_COLUMN_ORDER_SC2 =["Sample Name (from FASTA)", "Private ID", "GISAID ID (Public ID) - Optional", "Collection Date", "Collection Location", "Sequencing Date - Optional", "Sample is Private"]
    OUTPUT_COLUMN_ORDER_MPOX =["Sample Name (from FASTA)", "Private ID", "GenBank Accession (Public ID) - Optional", "Collection Date", "Collection Location", "Sequencing Date - Optional", "Sample is Private"]

    # set optional columns based on organism
    if "~{organism}" == "sars-cov-2":
      OPTIONAL_COLUMNS = OPTIONAL_COLUMNS_SC2
      OUTPUT_COLUMN_ORDER = OUTPUT_COLUMN_ORDER_SC2

    elif "~{organism}" == "mpox":
      OPTIONAL_COLUMNS = OPTIONAL_COLUMNS_MPOX
      OUTPUT_COLUMN_ORDER = OUTPUT_COLUMN_ORDER_MPOX

    # read in terra table and set the private_id_column_name to be only strings
    print("DEBUG: reading in terra table")
    table = pd.read_csv("~{full_terra_table}", delimiter='\t', header=0, dtype={"~{private_id_column_name}": 'str'})
    
    # extract the samples for upload from the entire table
    print("DEBUG: extracting samples from table")
    table = table[table["~{terra_table_name}_id"].isin("~{sep='*' sample_names}".split("*"))]

    # set all column headers to lowercase
    table.columns = table.columns.str.lower()

    # replace all NaN values with blanks
    table = table.fillna("")

    # extract metadata to new table
    print("DEBUG: copying require metadata columns")
    metadata = table[REQUIRED_COLUMNS].copy()

    print("DEBUG: copying optional metadata columns")
    for column in OPTIONAL_COLUMNS:
      if column in table.columns:
        metadata[column] = table[column]
      else:
        metadata[column] = ""
    
    metadata["Sample is Private"] = is_private

    metadata = metadata.fillna("")

    print("DEBUG: checking if private_id column was set")
    if "~{private_id_column_name}" == "~{terra_table_name}_id":
      print("DEBUG: removing duplicated column")
      metadata = metadata.loc[:, ~metadata.columns.duplicated()].copy()
      metadata["Private ID"] = metadata.loc[:, "~{private_id_column_name}"]
    else:
      metadata.rename(columns={"~{private_id_column_name}": "Private ID"}, inplace=True)

    # combine location data into one column
    print("DEBUG: combining collection location")
    metadata["Collection Location"] = metadata["~{continent_column_name}"] + "/" + metadata["~{country_column_name}"] + "/" + metadata["~{state_column_name}"]
    
    # add county to the location if the length of the county is > 0
    metadata["Collection Location"] = metadata.apply(lambda x: x["Collection Location"] + "/" + x["county"] if len(x["county"]) > 0 else x["Collection Location"], axis=1)

    if "~{organism}" == "sars-cov-2":
      print("DEBUG: creating GISAID ID column")
      metadata["GISAID ID (Public ID) - Optional"] = "hCoV-19/" + metadata["state"] + "/" + metadata["~{private_id_column_name}"] + "/" + metadata["collection_date"].apply(lambda x: year_getter(x))
    else:
      metadata.rename(columns={"~{genbank_accession_column_name}": "GenBank Accession (Public ID) - Optional"}, inplace=True)
    
    print("DEBUG: renaming the rest of the headers")
    # rename headers to match CZGenEpi's expected format
    metadata.rename(columns={"~{terra_table_name}_id": "Sample Name (from FASTA)", 
                             "~{collection_date_column_name}": "Collection Date",
                             "~{sequencing_date_column_name}": "Sequencing Date - Optional"}, inplace=True)

    # write the output to a csv file
    metadata.to_csv("czgenepi_prep_metadata.csv", columns=OUTPUT_COLUMN_ORDER, index=False)
    
    # create file transfer command and write it to a file
    table["cp"] = "gcloud storage cp " + table["~{assembly_fasta_column_name}"] + " ."
    table.to_csv("file-transfer.sh", columns=["cp"], index=False, header=False)

    # create fasta header renaming command and write it to a file
    table["rename_fasta_header"] = "sed -i '1s|.*|>" + table["~{terra_table_name}_id"] + "|' " + table.apply(lambda x: os.path.basename(x["~{assembly_fasta_column_name}"]), axis=1)
    table.to_csv("rename-command.sh", columns=["rename_fasta_header"], index=False, header=False)

    CODE

    echo "DEBUG: transfering fasta files"
    bash file-transfer.sh

    echo "DEBUG: renaming fasta headers"
    bash rename-command.sh

    echo "DEBUG: concatenating fasta files"
    cat *fasta > czgenepi_prep_concatenated.fasta

  >>>
  output {
    File concatenated_fasta = "czgenepi_prep_concatenated.fasta"
    File concatenated_metadata = "czgenepi_prep_metadata.csv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-08-08-2"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    #maxRetries: 3
  }
}
