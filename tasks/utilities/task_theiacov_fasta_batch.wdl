version 1.1

task sm_theiacov_fasta_wrangling { # the sm stands for supermassive
  input {
    String table_name
    String workspace_name
    String project_name
    String bucket_name

    Array[String] samplenames
    Map[String, File] sample_to_fasta
    String organism = "sars-cov-2"

    File? nextclade_tsv
    File? nextclade_json
    String? nextclade_docker
    String? nextclade_version
    String? nextclade_ds_tag

    File? pango_lineage_report
    String? pangolin_docker

    String seq_platform
    String assembly_method
    String theiacov_fasta_analysis_date
    String theiacov_fasta_version
    
    Int disk_size = 100
  }
  command <<<
    # convert the map into a JSON file for use in the python block
    # example map: {ERR4439752.test: /mnt/miniwdl_task_container/work/_miniwdl_inputs/0/ERR4439752.ivar.consensus.fasta}
    cp -v ~{write_json(sample_to_fasta)} sample_to_fasta.json
    
    # using "tmp" and will want to change that
    # this line splits into individual json files
    jq -c '.results = (.results[] | [.]) ' ~{nextclade_json} | awk '{ print > "tmp" NR ".json"}'
    # to-do: use gcloud storage cp to transfer to ~{bucket_name} location
    # split pangolin into individual csv files
    # and do the same thing with cp-ing
    # keep those filepaths, and save to table

    echo "DEBUG: Now entering Python block to perform parsing of data for populating sample-level table"

    python3 <<CODE 
    import pandas as pd 
    import numpy as np 
    import json
    import csv
    import os 
    import re
    
    # parse the map of sample names to fasta files
    with open("sample_to_fasta.json") as map_file:
      sample_to_fasta = json.load(map_file)
      # fix assembly_name
      print("trying to fix assembly names")
      sample_to_assembly = {name:re.split("[.]", os.path.basename(assembly))[0] for name, assembly in sample_to_fasta.items()}
      
    sample_name_array = "~{sep('*', samplenames)}".split("*")
    print(sample_name_array)


    # create a sample-level table to upload to terra
    upload_table = pd.DataFrame(sample_name_array, columns=["entity:~{table_name}_id"]).set_index("entity:~{table_name}_id")
    print(upload_table)

    upload_table["seq_platform"] = "~{seq_platform}"
    upload_table["assembly_method"] = "~{assembly_method}"
    upload_table["theiacov_fasta_analysis_date"] = "~{theiacov_fasta_analysis_date}"
    upload_table["theiacov_fasta_version"] = "~{theiacov_fasta_version}"

    # parse the NextClade output into an individual dataframe if a NextClade file exists
    if os.path.exists("~{nextclade_tsv}"):
      print("DEBUG: NEXTCLADE output TSV file identified; now parsing into a dataframe")
      nextclade = pd.read_csv("~{nextclade_tsv}", delimiter='\t')

      upload_table["nextclade_version"] = "~{nextclade_version}"
      upload_table["nextclade_docker"] = "~{nextclade_docker}"
      upload_table["nextclade_ds_tag"] = "~{nextclade_ds_tag}"
      
      for sample_name in sample_name_array:        
        print("DEBUG: the sample_name is {}".format(sample_name))
        assembly_name = sample_to_assembly[sample_name]
        print("DEBUG: the assembly_name is {}".format(assembly_name))

        if nextclade["seqName"].str.contains(assembly_name).any():
          if ("~{organism}" == "sars-cov-2"):
            nc_clade = str(nextclade.loc[nextclade["seqName"] == assembly_name]["clade_nextstrain"].item())
            who_clade = str(nextclade.loc[nextclade["seqName"] == assembly_name]["clade_who"].item())
            if (nc_clade != who_clade) and (nc_clade != "") and (who_clade != "") and (who_clade != "nan"):
              nc_clade = nc_clade + " (" + who_clade + ")"
            if nc_clade == "":
              nc_clade = "NA"
          else:
            nc_clade = str(nextclade.loc[nextclade["seqName"] == assembly_name]["clade"].item())
            if nc_clade == "":
              nc_clade = "NA"
          # replace nextclade value in datatable if exists, if not, create it
          if "nextclade_clade" not in upload_table.columns:
            upload_table["nextclade_clade"] = ""
          upload_table.at[sample_name, "nextclade_clade"] = nc_clade

          # parse nextclade_aa_subs
          nc_aa_subs = str(nextclade.loc[nextclade["seqName"] == assembly_name]["aaSubstitutions"].item())
          if nc_aa_subs == "":
            nc_aa_subs = "NA"
          elif ("~{organism}" == "flu"):
            print("FLU NOT SUPPORTED YET")
          if "nextclade_aa_subs" not in upload_table.columns:
            upload_table["nextclade_aa_subs"] = ""
          upload_table.at[sample_name, "nextclade_aa_subs"] = nc_aa_subs

          # parse nextclade_aa_dels
          nc_aa_dels = str(nextclade.loc[nextclade["seqName"] == assembly_name]["aaDeletions"].item())
          if nc_aa_dels == "":
            nc_aa_dels = "NA"
          if "nextclade_aa_dels" not in upload_table.columns:
            upload_table["nextclade_aa_dels"] = ""
          upload_table.at[sample_name, "nextclade_aa_dels"] = nc_aa_dels

          # parse nextclade_lineage
          try:
            nc_lineage = str(nextclade.loc[nextclade["seqName"] == assembly_name]["lineage"].item())
          except KeyError:
            nc_lineage = ""
          if nc_lineage == "":
            nc_lineage = "NA"
          if "nextclade_lineage" not in upload_table.columns:
            upload_table["nextclade_lineage"] = ""
          upload_table.at[sample_name, "nextclade_lineage"] = nc_lineage

    # parse the Pangolin lineage report into an individual dataframe if a Pangolin report file exists
    if os.path.exists("~{pango_lineage_report}"):
      print("DEBUG: PANGOLIN lineage report file identified; now parsing into a dataframe")
      pango_lineage_report = pd.read_csv("~{pango_lineage_report}", delimiter=',')
      
      upload_table["pangolin_docker"] = "~{pangolin_docker}"

      pangolin_version = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["pangolin_version"].item()
      version = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["version"].item()
      upload_table["panglin_version"] = "pangolin {}; {}".format(pangolin_version, version)

      # iterate through results and add to table
      for sample_name in sample_name_array:        
        assembly_name = sample_to_assembly[sample_name]
 
        # parse pango_lineage from pango lineage report
        pango_lineage = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["lineage"].item()
        if "pango_lineage" not in upload_table.columns:
          upload_table["pango_lineage"] = ""
        upload_table.at[sample_name, "pango_lineage"] = pango_lineage

        # parse pango_lineage_expanded from pango lineage report
        try:
          pango_lineage_expanded = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["expanded_lineage"].item()
        except KeyError:
            pango_lineage_expanded = ""
        if "pango_lineage_expanded" not in upload_table.columns:
          upload_table["pango_lineage_expanded"] = ""
        upload_table.at[sample_name, "pango_lineage_expanded"] = pango_lineage_expanded

        # parse pangolin_conflicts from pango lineage report
        pangolin_conflicts = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["conflict"].item()
        if "pangolin_conflicts" not in upload_table.columns:
          upload_table["pangolin_conflicts"] = ""
        upload_table.at[sample_name, "pangolin_conflicts"] = pangolin_conflicts

        # parse pangolin_notes from pango lineage report
        pangolin_notes = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["note"].item()
        if "pangolin_notes" not in upload_table.columns:
          upload_table["pangolin_notes"] = ""
        upload_table.at[sample_name, "pangolin_notes"] = pangolin_notes
    
    # TODO: drop assembly_fasta column? What about other columns that are not required for upload?
    upload_table.to_csv("TERRA_TABLE_TEMP.tsv", sep='\t', index=True)

    CODE



    # upload results to terra databable 
    #python3 scripts/import_large_tsv/import_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --tsv <path_to_tsv_to_upload> 

    echo "DEBUG: upload to terra table complete"
  >>>
  output {
    File terra_table = "TERRA_TABLE_TEMP.tsv"
    Boolean success = true
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 4
    disks:  "~{disk_size} local-disk "
    disk: "~{disk_size} GB"
    preemptible: 0
  }
}