version 1.0

task terra_to_bigquery {
  input {
    Array[String]  terra_projects
    Array[String]  workspace_names
    Array[String]  table_names
    Array[String]  table_ids
    Array[String]  gcs_uri_prefixs
    Array[String]  output_filename_prefix
    String  docker = "us-docker.pkg.dev/general-theiagen/broadinstitute/terra-tools:tqdm"
    Int page_size = 5000
    Int mem_size_gb = 32
    Int CPUs = 8
    Int disk_size = 100
    String sleep_time = "15m"
  }
  meta {
    volatile: true
  }
  command <<<
  set -e

  # set bash arrays
  terra_project_array=(~{sep=' ' terra_projects})
  terra_project_array_len=$(echo "${#terra_project_array[@]}")
  workspace_name_array=(~{sep=' ' workspace_names})
  workspace_name_array_len=$(echo "${#workspace_name_array[@]}")
  table_name_array=(~{sep=' ' table_names})
  table_name_array_len=$(echo "${#table_name_array[@]}")
  table_id_array=(~{sep=' ' table_ids})
  table_id_array_len=$(echo "${#table_id_array[@]}")
  gcs_uri_prefix_array=(~{sep=' ' gcs_uri_prefixs})
  gcs_uri_prefix_array_len=$(echo "${#gcs_uri_prefix_array[@]}")
  output_filename_prefix_array=(~{sep=' ' output_filename_prefix})
  output_filename_prefix_array_len=$(echo "${#output_filename_prefix_array[@]}")

  # Ensure equal length of all input arrays
  echo "Terra Projects array length: $terra_project_array_len"
  echo "Workspace name array length: $workspace_name_array_len"
  echo "Table Name array length: $table_name_array_len"
  echo "Table ID array length: $table_id_array_len" 
  echo "GCS URI prefixes array length: $gcs_uri_prefix_array_len"
  echo "Output filename prefix array: $output_filename_prefix_array_len"
  echo

  echo "comparing arrays: terra_project_array_len and workspace_name_array_len"
  if [ $terra_project_array_len -eq $workspace_name_array_len ]; then
    echo "Input arrays are of equal length."
  else 
    echo "Input arrays are of unequal length. Exiting"
    exit 1
  fi 
  
  echo "comparing arrays: terra_project_array_len and table_name_array_len"
  if [ $terra_project_array_len -eq $table_name_array_len ]; then
    echo "Input arrays are of equal length."
  else 
    echo "Input arrays are of unequal length. Exiting. Please check your inputs!"
    exit 1
  fi
  
  echo "comparing arrays: terra_project_array_len and table_id_array_len"
  if [ $terra_project_array_len -eq $table_id_array_len ]; then
    echo "Input arrays are of equal length."
  else 
    echo "Input arrays are of unequal length. Exiting. Please check your inputs!"
    exit 1
  fi
  
  echo "comparing arrays: terra_project_array_len and gcs_uri_prefix_array_len"
  if [ $terra_project_array_len -eq $gcs_uri_prefix_array_len ]; then
    echo "Input arrays are of equal length."
  else 
    echo "Input arrays are of unequal length. Exiting. Please check your inputs!"
    exit 1
  fi
  
   echo "comparing arrays: terra_project_array_len and output_filename_prefix_array_len"
  if [ $terra_project_array_len -eq $output_filename_prefix_array_len ]; then
    echo "Input arrays are of equal length."
  else 
    echo "Input arrays are of unequal length. Exiting. Please check your inputs!"
    exit 1
  fi
  
  # [ $terra_project_array_len -eq $output_filename_prefix_array_len ]; then
  #   echo "Input arrays are of equal length."
  #   echo "Proceeding to transfer the following Terra Data Tables to their specified GCS URIs: ${#gcs_uri_prefix_array[@]}"
  #   echo
  #   echo "Table name array: ${#table_name_array[@]}"
  #   echo
  #   echo "Transfer will occur every ~{sleep_time} until this job is aborted."
  # else
  #   echo "Input arrays are of unequal length" >&2
  #   echo "Terra Projects array length: $terra_project_array_len" >&2
  #   echo "Workspace name array length: $workspace_name_array_len"  >&2
  #   echo "Table Name array length: $table_name_array_len" >&2
  #   echo "Table ID array length: $table_id_array_len" >&2
  #   echo "GCS URI prefix array length: $gcs_uri_prefix_array_len" >&2
  #   echo "Output filename prefix array: $output_filename_prefix_array_len" >&2
  #   echo
  #   echo "Exiting script! Please check your inputs!"
  #   exit 1
  # fi

  # Infinite While loop
  counter=0
  echo -e "**ENTERING LOOP**"
  while true
  do

    # counter and sanity checks for troubleshooting
    counter=$((counter+1))
    date_tag=$(date +"%Y-%m-%d-%Hh-%Mm-%Ss")
    echo -e "\n========== Iteration number ${counter} of continuous loop =========="
    echo "TIME: ${date_tag}"

  # Loop through inputs and run python script to create tsv/json and push json to specified gcp bucket
    for index in "${!terra_project_array[@]}"; do
      date_tag=$(date +"%Y-%m-%d-%Hh-%Mm-%Ss")
      terra_project=${terra_project_array[$index]}
      workspace_name=${workspace_name_array[$index]}
      table_name=${table_name_array[$index]}
      # re-enabling table_id to use to fill out the "source_terra_table" field in output JSON instead of using the table_name 
      # allows to differentiate between identically-named data tables coming from different counties or labs. eg. "sample" data table for all ClearLabs users
      table_id=${table_id_array[$index]}
      gcs_uri=${gcs_uri_prefix_array[$index]}
      output_filename_prefix=${output_filename_prefix_array[$index]}

      # if user specifies 'date' as output_filename_prefix, reset variable to default table filename ${table_name}_${date_tag}.json
      if [ "${output_filename_prefix}" == "date" ]; then
        echo 'User specified "date"' "for output_filename_prefix, final output JSON will be named: ${gcs_uri}${table_name}_${date_tag}.json"
        output_filename_prefix="${table_name}_${date_tag}"
      fi

      export terra_project workspace_name table_name table_id date_tag gcs_uri output_filename_prefix

      echo
      echo "***Exporting Terra table ${table_name} from workspace ${workspace_name} in Terra project ${terra_project}***"

      # download Terra table TSV using export_large_tsv.py from Broad
      python3 /scripts/export_large_tsv/export_large_tsv.py \
        --project "${terra_project}" \
        --workspace "${workspace_name}" \
        --entity_type "${table_name}" \
        --page_size ~{page_size} \
        --tsv_filename "${table_name}_${date_tag}.tsv"

      echo -e "\n::Procesing ${table_name} for export (${date_tag})::"
      echo
      echo "entering python block of code...."

      # reformat TSV using code below
      # additionally take cleaned-TSV and create nlJSON
      python3<<CODE
  import csv
  import json
  import collections
  import os

  from firecloud import api as fapi

  # sanity checks for env variables loaded into python
  workspace_project = os.environ['terra_project']
  print("workspace project: "+ workspace_project)
  workspace_name = os.environ['workspace_name']
  print("workspace name: "+ workspace_name)
  table_name = os.environ['table_name']
  print("table name: "+ table_name)
  table_id = os.environ['table_id']
  print("table_id: "+ table_id)
  #out_fname = os.environ['table_id']
  #print("out_fname: " + out_fname)
  date_tag = os.environ['date_tag']
  print("date_tag: " + date_tag)

  ####COMMENTING OUT THIS BLOCK####
  # Grabbbing defined table using firecloud api and reading data to to python dictionary
  #table = json.loads(fapi.get_entitiesget_entities(workspace_project, workspace_name, table_name).text)

  # instead of loading JSON directly from Terra data table, load in TSV that was just exported from Terra

  ####COMMENTING OUT THIS BLOCK####
  # This block transforms JSON to dictionary
  # headers = collections.OrderedDict()
  # rows = []
  # headers[table_name + "_id"] = 0
  # for row in table:
  #   outrow = row['attributes']
  #   for x in outrow.keys():
  #     headers[x] = 0
  #     if type(outrow[x]) == dict and set(outrow[x].keys()) == set(('itemsType', 'items')):
  #       outrow[x] = outrow[x]['items']
  #   outrow[table_name + "_id"] = row['name']
  #   rows.append(outrow)

  ####COMMENTING OUT THIS BLOCK####
  # Writing tsv output from dictionary object
  # with open(out_fname+'_temp.tsv', 'w') as outf:
  #   writer = csv.DictWriter(outf, headers.keys(), delimiter='\t', dialect=csv.unix_dialect, quoting=csv.QUOTE_MINIMAL)
  #   writer.writeheader()
  #   writer.writerows(rows)
  
  print("adding source_terra_table column to TSV...")

  # TSV add additional column
  # Add column to capture source terra table (table_name) 
  with open(table_name + '_' + date_tag +'.tsv','r') as csvinput:
    with open(table_name+'.tsv', 'w') as csvoutput:
        writer = csv.writer(csvoutput, delimiter='\t')
        reader = csv.reader(csvinput, delimiter='\t')

        all = []
        tsv_row = next(reader)
        tsv_row.append("source_terra_table")
        all.append(tsv_row)

        for tsv_row in reader:
            tsv_row.append(table_id)
            all.append(tsv_row)

        writer.writerows(all)

  print("converting TSV to newline JSON...")

  # Writing the newline json file from tsv output above
  with open(table_name+'.tsv', 'r') as infile:
    headers = infile.readline()
    headers_array = headers.strip().split('\t')
    headers_array[0] = "specimen_id"
    with open(table_name+'.json', 'w') as outfile:
      for line in infile:
        outfile.write('{')
        line_array=line.strip().split('\t')
        for x,y in zip(headers_array, line_array):
          if x == "nextclade_aa_dels" or x == "nextclade_aa_subs":
            y = y.replace("|", ",")
          if y == "NA":
            y = ""
          if y == "N/A":
            y = ""
          if y == "Unknown":
            y = ""
          if y == "unknown":
            y = ""
          if y == "UNKNOWN":
            y = ""
          if y == "required_for_submission":
            y = ""
          if "Uneven pairs:" in y:
            y = ""
          if x == "County":
            pass
          else:
            outfile.write('"'+x+'"'+':'+'"'+y+'"'+',')
        outfile.write('"notes":""}'+'\n')
  print ("finished creating newline JSON, exiting python block...")
  CODE

      export CLOUDSDK_PYTHON=python2.7  # ensure python 2.7 for gsutil commands

      # add date tag when transferring file to gcp
      #### date_tag variable is already set above the python block, so commenting out ###
      #date_tag=$(date +"%Y-%m-%d-%Hh-%Mm-%Ss")

      # THIS IF BLOCK SHOULD ALWAYS BE TRIGGERED AS LONG AS USER DEFINES OUTPUT FILENAME PREFIX FOR ALL TABLES
      # if user defines a filename prefix, then use it to name the output JSON file
      # if output_filename_prefix bash input string is non-zero, return TRUE
      if [ -n "${output_filename_prefix}" ]; then
        echo
        echo "User specified an output filename prefix of: ${output_filename_prefix}"
        # copy new line JSON to bucket & copy re-formatted TSV (for testing purposes)
        gsutil -m cp "${table_name}.json" "${gcs_uri}${output_filename_prefix}.json"
        echo "${output_filename_prefix}.json copied to ${gcs_uri}"
      else
        # copy new line JSON to bucket & copy re-formatted TSV (for testing purposes)
        echo "User did NOT specify an output prefix, using default prefix with table_name and date_tag variables"
        echo
        gsutil -m cp "${table_name}.json" "${gcs_uri}${table_name}_${date_tag}.json"
        echo "${table_name}_${date_tag}.json copied to ${gcs_uri}"
        echo
      fi
      
      echo -e "***Finished exporting and copying of table ${table_name} to the specified google bucket. Moving to the next data table!***"
      echo

      unset CLOUDSDK_PYTHON   # probably not necessary, but in case I do more things afterwards, this resets that env var
    done
    echo "Sleeping for user-specified time of " ~{sleep_time}
    sleep ~{sleep_time}
    echo "Finished sleeping, onto the next iteration of the loop!"
  done
  echo "Loop exited"
  >>>
  runtime {
    docker: docker
    memory: "~{mem_size_gb} GB"
    cpu: CPUs
    disks: "local-disk ~{disk_size} SSD"
  }
  output {
    ## add outputs for all intermediate files
  }
}