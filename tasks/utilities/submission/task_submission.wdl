version 1.0

task prune_table {
  input {
    String table_name
    String workspace_name
    String project_name
    File? input_table
    Array[String] sample_names
    String biosample_type
    String bioproject
    String gcp_bucket_uri
    Boolean skip_biosample
    String read1_column_name = "read1"
    String read2_column_name = "read2"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --entity_type ~{table_name} --tsv_filename ~{table_name}-data.tsv
    
    # when running locally, use the input_table in place of downloading from Terra
    #cp ~{input_table} ~{table_name}-data.tsv

    # transform boolean skip_biosample into string for python comparison
    if ~{skip_biosample}; then
      export skip_bio="true"
    else 
      export skip_bio="false"
    fi

    python3 <<CODE 
    import pandas as pd
    import numpy as np
    import os

    # set a function to remove NA values and return the cleaned table and a table of excluded samples
    def remove_nas(table, required_metadata):
      table.replace(r'^\s+$', np.nan, regex=True) # replace blank cells with NaNs 
      excluded_samples = table[table[required_metadata].isna().any(axis=1)] # write out all rows that are required with NaNs to a new table
      excluded_samples.set_index("~{table_name}_id".lower(), inplace=True) # convert the sample names to the index so we can determine what samples are missing what
      excluded_samples = excluded_samples[excluded_samples.columns.intersection(required_metadata)] # remove all optional columns so only required columns are shown
      excluded_samples = excluded_samples.loc[:, excluded_samples.isna().any()] # remove all NON-NA columns so only columns with NAs remain; Shelly is a wizard and I love her 
      table.dropna(subset=required_metadata, axis=0, how='any', inplace=True) # remove all rows that are required with NaNs from table

      return table, excluded_samples

    # read export table into pandas
    tablename = "~{table_name}-data.tsv"
    table = pd.read_csv(tablename, delimiter='\t', header=0, dtype={"~{table_name}_id": 'str'}) # ensure sample_id is always a string)

    # extract the samples for upload from the entire table
    table = table[table["~{table_name}_id"].isin("~{sep='*' sample_names}".split("*"))]

    # set required and optional metadata fields based on the biosample_type package
    if ("~{biosample_type}".lower() == "microbe"):
      required_metadata = ["submission_id", "organism", "collection_date", "geo_loc_name", "sample_type"]
      optional_metadata = ["sample_title", "bioproject_accession", "attribute_package", "strain", "isolate", "host", "isolation_source", "altitude", "biomaterial_provider", "collected_by", "depth", "env_broad_scale", "genotype", "host_tissue_sampled", "identified_by", "lab_host", "lat_lon", "mating_type", "passage_history", "samp_size", "serotype", "serovar", "specimen_voucher", "temp", "description", "MLST"]
      # add a column for biosample package -- required for XML submission
      table["attribute_package"] = "Microbe.1.0"
      # future qc checks:
      #   q-score >= 30
      #   reads > 50 bp
      #   trailing/leading bases removed
      #   similar GC content to expected genome
      #   assembled genome ratio ~1.0
      #   200 contigs or less

    elif ("~{biosample_type}".lower() == "wastewater"):
      required_metadata = ["submission_id", "organism", "collection_date", "geo_loc_name", "isolation_source", "ww_population", "ww_sample_duration", "ww_sample_matrix", "ww_sample_type", "ww_surv_target_1", "ww_surv_target_1_known_presence"]
      optional_metadata = ["sample_title", "bioproject_accession", "attribute_package", "collected_by", "purpose_of_ww_sampling","purpose_of_ww_sequencing", "sequenced_by", "ww_endog_control_1", "ww_endog_control_1_conc", "ww_endog_control_1_protocol", "ww_endog_control_1_units", "ww_endog_control_2", "ww_endog_control_2_conc", "ww_endog_control_2_protocol", "ww_endog_control_2_units", "ww_flow", "ww_industrial_effluent_percent", "ww_ph", "ww_population_source", "ww_pre_treatment", "ww_primary_sludge_retention_time", "ww_processing_protocol", "ww_sample_salinity", "ww_sample_site", "ww_surv_jurisdiction", "ww_surv_system_sample_id", "ww_surv_target_1_conc", "ww_surv_target_1_conc_unit", "ww_surv_target_1_extract", "ww_surv_target_1_extract_unit", "ww_surv_target_1_gene", "ww_surv_target_1_protocol", "ww_surv_target_2", "ww_surv_target_2_conc", "ww_surv_target_2_conc_unit", "ww_surv_target_2_extract", "ww_surv_target_2_extract_unit", "ww_surv_target_2_gene", "ww_surv_target_2_known_present", "ww_surv_target_2_protocol", "ww_temperature", "ww_total_suspended_solids", "description"]
      # add a column for biosample package -- required for XML submission
      table["attribute_package"] = "SARS-CoV-2.wwsurv.1.0"

      # qc checks:

   
    elif ("~{biosample_type}".lower() == "pathogen") or ("pathogen.cl" in "~{biosample_type}".lower()):
      required_metadata = ["submission_id", "organism", "collected_by", "collection_date", "geo_loc_name", "host", "host_disease", "isolation_source", "lat_lon"]
      optional_metadata = ["sample_title", "isolation_type", "bioproject_accession", "attribute_package", "strain", "isolate", "culture_collection", "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state", "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar", "specimen_voucher", "subgroup", "subtype", "description"] 
      # add a column for biosample package -- required for XML submission
      table["attribute_package"] = "Pathogen.cl"
      # future qc checks:
      #   gc after trimming 42-47.5%
      #   average phred after trimming >= 28
     
      #   coverage after trimming >= 20X
      #if "mean_coverage_depth" in table.columns:
      #  table = table[(table.mean_coverage_depth > 20)]
    elif ("pathogen.env" in "~{biosample_type}".lower()):
      required_metadata = ["submission_id", "organism", "collected_by", "collection_date", "geo_loc_name", "isolation_source", "lat_lon"]
      optional_metadata = ["host", "host_disease", "isolation_type", "sample_title", "bioproject_accession", "attribute_package", "strain", "isolate", "culture_collection", "genotype", "host_age", "host_description", "host_disease_outcome", "host_disease_stage", "host_health_state", "host_sex", "host_subject_id", "host_tissue_sampled", "passage_history", "pathotype", "serotype", "serovar", "specimen_voucher", "subgroup", "subtype", "description"] 
      # add a column for biosample package -- required for XML submission
      table["attribute_package"] = "Pathogen.env.1.0"
      # future qc checks:
      #   gc after trimming 42-47.5%
      #   average phred after trimming >= 28
     
      #   coverage after trimming >= 20X
      #if "mean_coverage_depth" in table.columns:
      #  table = table[(table.mean_coverage_depth > 20)]

    elif ("~{biosample_type}".lower() == "virus"):
      required_metadata = ["submission_id", "organism", "isolate", "collection_date", "geo_loc_name", "isolation_source"]
      optional_metadata = ["sample_title", "bioprojection_accession","attribute_package", "host", "lab_host", "altitude", "biomaterial_provider", "collected_by", "culture_collection", "depth", "disease", "env_broad_scale", "genotype", "host_tissue_sampled", "identified_by", "lat_lon", "passage_history", "samp_size", "serotype", "specimen_voucher", "strain", "temp", "description"]

      table["attribute_package"] = "Virus.1.0"


    else:
      raise Exception('Only "Microbe", "Virus", "Pathogen" and "Wastewater" are supported as acceptable input for the \`biosample_type\` variable at this time. You entered ~{biosample_type}.')

    # sra metadata is the same regardless of biosample_type package, but I'm separating it out in case we find out this is incorrect
    sra_required = ["~{table_name}_id", "submission_id", "library_ID", "title", "library_strategy", "library_source", "library_selection", "library_layout", "platform", "instrument_model", "design_description", "filetype", "~{read1_column_name}"]
    sra_optional = ["~{read2_column_name}"]

    # if biosample accessions are provided, add those to the end of the sra_required field
    if (os.environ["skip_bio"] == "true"):
      sra_required.append("biosample_accession")

    # combine all required fields into one array for easy removal of NaN cells
    required_fields = required_metadata + sra_required

    # remove required rows with blank cells from table
    table, excluded_samples = remove_nas(table, required_fields)
    with open("excluded_samples.tsv", "a") as exclusions:
      exclusions.write("Samples excluded for missing required metadata (will have empty values in indicated columns):\n")
    excluded_samples.to_csv("excluded_samples.tsv", mode='a', sep='\t')

    # add bioproject_accesion to table
    table["bioproject_accession"] = "~{bioproject}"
    
    # extract the required metadata from the table
    biosample_metadata = table[required_metadata].copy()

    # add optional metadata fields if present; rename first column
    for column in optional_metadata:
      if column in table.columns:
        biosample_metadata[column] = table[column]
    biosample_metadata.rename(columns={"submission_id" : "sample_name"}, inplace=True)

    # extract the required metadata from the table; rename first column 
    sra_metadata = table[sra_required].copy()
    for column in sra_optional:
      if column in table.columns:
        sra_metadata[column] = table[column]
    sra_metadata.rename(columns={"submission_id" : "sample_name"}, inplace=True)

    # prettify the filenames and rename them to be sra compatible
    sra_metadata["~{read1_column_name}"] = sra_metadata["~{read1_column_name}"].map(lambda filename: filename.split('/').pop())
    sra_metadata.rename(columns={"~{read1_column_name}" : "filename"}, inplace=True)
    table["~{read1_column_name}"].to_csv("filepaths.tsv", index=False, header=False) # make a file that contains the names of all the reads so we can use gsutil -m cp
    if "~{read2_column_name}" in sra_metadata.columns:
      sra_metadata["~{read2_column_name}"] = sra_metadata["~{read2_column_name}"].map(lambda filename2: filename2.split('/').pop())   
      sra_metadata.rename(columns={"~{read2_column_name}" : "filename2"}, inplace=True)
      table["~{read2_column_name}"].to_csv("filepaths.tsv", mode='a', index=False, header=False)
    
    # write metadata tables to tsv output files
    biosample_metadata.to_csv("biosample_table.tsv", sep='\t', index=False)
    sra_metadata.to_csv("sra_table_to_edit.tsv", sep='\t', index=False)

    CODE

    # prune the first two columns of sra_table_to_edit to remove the tablename_id and submission_id columns
    cut -f3- sra_table_to_edit.tsv > sra_table.tsv

    # copy the raw reads to the bucket specified by user
    export CLOUDSDK_PYTHON=python2.7  # ensure python 2.7 for gsutil commands
    # iterate through file created earlier to grab the uri for each read file
    while read -r line; do
      echo "running \`gsutil -m cp ${line} ~{gcp_bucket_uri}\`"
      gsutil -m cp -n ${line} ~{gcp_bucket_uri}
    done < filepaths.tsv
    unset CLOUDSDK_PYTHON   # probably not necessary, but in case I do more things afterwards, this resets that env var

  >>>
  output {
    File biosample_table = "biosample_table.tsv"
    File sra_table = "sra_table.tsv"
    File sra_table_for_biosample = "sra_table_to_edit.tsv"
    File excluded_samples = "excluded_samples.tsv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task add_biosample_accessions {
  input {
    File generated_accessions
    File sra_metadata
    String project_name
    String workspace_name
    String table_name
  }
  command <<<
    echo "Uploading biosample_accession to the Terra data table"

    ## check if any biosample accessions were made 
    tail -n +2 ~{generated_accessions} > removed_header
    if [ -s removed_header ]; then
      echo true > PROCEED
      echo "Biosample accessions were generated! Proceeding to SRA submission" > BIOSAMPLE_STATUS
      
      # add biosample accessions to sra_metadata table:

      # extract the table_id column from sra_metadata and the biosample accession from attributes and output to table
      awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} {print $1, a[$2]}' ~{generated_accessions} ~{sra_metadata} > table-ids-and-biosamples.tsv
      #  BEGIN {OFS="\t"} sets the output field separator to tab
      #  FNR==NR is an if statement saying that you do the first command when true, and the second when false
      #     this is false when attributes is finished being read and the sra_metadata file starts being read
      #  {a[$2]=$1; next} creates and array where the 2nd column of generated_accessions (sample_name) is the index 
      #     and is set equal to column 1 (biosample accession)
      #  {print $1, a[$2]} prints the table_id column and then the biosample_accession in the array that matches
      #     the second column of sra_metadata (sample_name)

      # echo out the header for the upload table
      echo -e "entity:~{table_name}_id\tbiosample_accession" > upload-terra.tsv

      # skip the header and append to new file
      tail -n +2 table-ids-and-biosamples.tsv >> upload-terra.tsv

      # upload biosample_accessions to Terra
      python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --tsv upload-terra.tsv

      echo "Adding biosample_accession to the sra_metadata table"
      # extract from the attributes file the biosample and original name columns
      # put the original name in column 1, biosample in column 2
      awk -F '\t' '{print $2, $1}' OFS='\t' ~{generated_accessions} > biosample_temp.tsv

      # remove the table_id column
      cut -f2- ~{sra_metadata} > sra_temp.tsv

      # echo out the header for the updated sra_metadata file
      echo -e "$(head -n 1 sra_temp.tsv)\tbiosample_accession" > "sra_table_with_biosample_accessions-with-sample-names.tsv"

      # join the biosample_temp with the sra_metadata; using tail to skip the header 
      # this join will remove any sra numbers do not have biosample accessions
      join -t $'\t' <(sort <(tail -n+2 sra_temp.tsv)) <(sort <(tail -n+2 biosample_temp.tsv)) >> "sra_table_with_biosample_accessions-with-sample-names.tsv"

      # remove the unnecessary submission_id column 
      cut -f2- "sra_table_with_biosample_accessions-with-sample-names.tsv" > "sra_table_with_biosample_accessions.tsv"

    else # no biosample_accessions generated
      echo false > PROCEED
      echo "Biosample accessions were NOT generated! Ending process now" > BIOSAMPLE_STATUS
    fi

  >>>
  output {
    File? sra_table = "sra_table_with_biosample_accessions.tsv"
    String biosample_status = read_string("BIOSAMPLE_STATUS")
    Boolean proceed = read_boolean("PROCEED")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}