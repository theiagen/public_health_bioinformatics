version 1.0

task summarize_data {
  input {
    Array[String]? sample_names
    String? terra_project
    String? terra_workspace
    String? terra_table
    String? column_names # string of comma-delimited column names
    String? output_prefix
    String? id_column_name

    Int disk_size = 100
    Int cpu = 8
    Int memory = 1 
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16" 
    # commenting out this option since it's for local dev. Prefer this option to not appear in Terra
    #File? input_table
    Boolean phandango_coloring = false
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<   
    set -euo pipefail
    
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --entity_type ~{terra_table} --tsv_filename ~{terra_table}-data.tsv 
    
    # when running locally, use the input_table in place of downloading from Terra
    # TO RENABLE: uncomment line below, and add back tilde in front of {input_table}
    #cp {input_table} ~{terra_table}-data.tsv
    
    if ~{phandango_coloring}; then
      export phandango_coloring="true"
    else 
      export phandango_coloring="false"
    fi

    # indicate if a different id_column should be used than the default
    if [[ -z "~{id_column_name}" ]]; then
      export default_column="true"
    else
      export default_column="false"
    fi

    python3 <<CODE 
  import pandas as pd
  import numpy as np
  import itertools
  import os
  import re

  # read exported Terra table into pandas
  tablename = "~{terra_table}-data.tsv" 
  table = pd.read_csv(tablename, delimiter='\t', header=0, index_col=False, dtype={"~{terra_table}_id": 'str'}) # ensure sample_id is always a string

  # extract the samples for upload from the entire table
  if (os.environ["default_column"] == "true"):
    table = table[table["~{terra_table}_id"].isin("~{sep='*' sample_names}".split("*"))]
  else:
    table = table[table["~{id_column_name}"].isin("~{sep='*' sample_names}".split("*"))] 

  # cast entire table as str
  table = table.astype(str)
  
  # split comma-separated column list into an array
  columns = "~{column_names}".split(",")

  # incorporate ID column
  temporarylist = []
  if (os.environ["default_column"] == "true"):
    id_col = "~{terra_table}_id"
    temporarylist.append(id_col)
  else:
    id_col = "~{id_column_name}"
    temporarylist.append(id_col)
  temporarylist += columns

  # only represent extracted columns
  table = table[temporarylist].copy()

  # create a table to search through containing only columns of interest
  searchtable = table[columns].copy()

  # set index and output to tsv w/index
  if (os.environ["default_column"] == "true"):
    filteredmetadata = searchtable.set_index(table["~{terra_table}_id"])
  else:
    filteredmetadata = searchtable.set_index(table["~{id_column_name}"])
  filteredmetadata.to_csv("~{output_prefix}_filtered_metadata.tsv", sep='\t', index=True)

  # iterate through the columns of interest and combine into a single list
  genes = []
  for item in columns:
    genes.append(table[item].str.split(",").explode().tolist())

  # add phandango coloring tags if indicated
  if (os.environ["phandango_coloring"] == "true"):
    i = 1 # initialize a starting group at 1
    newgenes = [] # create a new temporary list
    for group in genes: # iterate through gene list by items found within a column
      if (i < 10):
        newgroup = [] # create a new temporary sublist
        for item in group: # for every item in a column, 
          newitem = str(item) + ":o" + str(i) # add a unique :o coloring (as indicated by the str(i))
          newgroup.append(newitem) # add phandango-suffixed item to the new sublist
        newgenes.append(newgroup) # add the new sublist to the new list
        i += 1 # increment the i value so each column gets its own coloring     
      else:
        print("DEBUG: Some items will be ignored because a maximum of ten columns can be presented at once with phandango coloring.")

    # overwrite genes with newgenes (which now has the phandango coloring suffix)
    genes = newgenes 
  else:
    print("DEBUG: Phandango coloring was not applied")

  # flattening the list
  genes = list(itertools.chain.from_iterable(genes))

  # removing duplicates but maintaining order
  genes = sorted(set(genes), key=lambda x: genes.index(x))

  # force string type for pattern matching
  searchtable = searchtable.astype(str)

  # add genes as true/false entries into table row-by-row
  # iterate from largest to smallest string size to prevent substrings from matching large strings
  for item in sorted(genes, key=lambda x: len(x), reverse=True):
    if (os.environ["phandango_coloring"] == "true"): # remove coloring suffix (CAUTION: ASSUMES LESS THAN 10 COLUMN NAMES PROVIDED)
      table[item] = searchtable.apply(lambda row: row.str.contains(re.escape(item[:len(item)-3])).any(), axis=1)
      # remove the entry to prevent future pattern matching
      searchtable = searchtable.apply(lambda row: row.str.replace(re.escape(item[:len(item)-3]), ""))
    else:
      table[item] = searchtable.apply(lambda row: row.str.contains(re.escape(item)).any(), axis=1)
      searchtable = searchtable.apply(lambda row: row.str.replace(re.escape(item), ""))

  # sort table by original gene index
  final_cols = [id_col] + genes
  final_table = table[final_cols]

  # replace all "False" cells with null values 
  final_table[final_table.eq(False)] = np.nan

  final_table.to_csv("~{output_prefix}_summarized_data.csv", sep=',', index=False)

  CODE
  >>>
  output {
    File summarized_data = "~{output_prefix}_summarized_data.csv"
    File filtered_metadata = "~{output_prefix}_filtered_metadata.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}