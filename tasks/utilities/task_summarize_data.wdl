version 1.0

task summarize_data {
  input {
    Array[String]? sample_names
    String? terra_project
    String? terra_workspace
    String? terra_table
    String? column_names # string of space delimited column names
    
    Int disk_size = 100
    File? input_table
    Boolean phandango_coloring = true
  }
  command <<<   
    # when running on terra, comment out all input_table mentions
    python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --entity_type ~{terra_table} --tsv_filename ~{terra_table}-data.tsv 
    
    # when running locally, use the input_table in place of downloading from Terra
    #cp ~{input_table} ~{terra_table}-data.tsv
    
    if ~{phandango_coloring}; then
      export phandango_coloring="true"
    else 
      export phandango_coloring="false"
    fi

    python3 <<CODE 
  import pandas as pd
  import itertools
  import os
  import re

  # read exported Terra table into pandas
  tablename = "~{terra_table}-data.tsv" 
  table = pd.read_csv(tablename, delimiter='\t', header=0, dtype={"~{terra_table}_id": 'str'}) # ensure sample_id is always a string

  # extract the samples for upload from the entire table
  table = table[table["~{terra_table}_id"].isin("~{sep='*' sample_names}".split("*"))]

  # split column list into an array
  columns = "~{column_names}".split(" ")

  temporarylist = []
  temporarylist.append("~{terra_table}_id")
  temporarylist += columns

  table = table[temporarylist].copy()

  # create a table to search through containing only columns of interest
  searchtable = table[columns].copy()

  # iterate through the columns of interest and combine into a single list
  genes = []
  for item in columns:
    genes.append(table[item].str.split(",").explode().tolist())

  # add phandango coloring tags if indicated
  if (os.environ["phandango_coloring"] == "true"):
    i = 1 # initialize a starting group at 1
    newgenes = [] # create a new temporary list
    for group in genes: # iterate through gene list by items found within a column
      newgroup = [] # create a new temporary sublist
      for item in group: # for every item in a column, 
        newitem = item + ":o" + str(i) # add a unique :o coloring (as indicated by the str(i))
        newgroup.append(newitem) # add phandango-suffixed item to the new sublist
      newgenes.append(newgroup) # add the new sublist to the new list
      i += 1 # increment the i value so each column gets its own coloring     

    # overwrite genes with newgenes (which now has the phandango coloring suffix)
    genes = newgenes 
  else:
    print("NOTE: Phandango coloring was not applied")

  # flattening the list
  genes = list(itertools.chain.from_iterable(genes))

  # removing duplicates
  genes = list(set(genes))

  # add genes as true/false entries into table
  for item in genes:
    if (os.environ["phandango_coloring"] == "true"): # remove coloring suffix (CAUTION: ASSUMES LESS THAN 10 COLUMN NAMES PROVIDED)
      table[item] = searchtable.apply(lambda row: row.astype(str).str.contains(re.escape(item[:len(item)-3])).any(), axis=1)
    else:
      table[item] = searchtable.apply(lambda row: row.astype(str).str.contains(re.escape(item)).any(), axis=1)

  # dropping columns of interest so only true/false ones remain
  table.drop(columns,axis=1,inplace=True)

  table.to_csv("summarized_data.csv", sep=',', index=False)

  CODE
  >>>
  output {
    File summarized_data = "summarized_data.csv"
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}