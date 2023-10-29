version 1.0

task sm_theiacov_fasta_wrangling { # the sm stands for supermassive
  input {
    String table_name
    String workspace_name
    String project_name
    File? input_table
    Array[String] samplenames
    File? nextclade_tsv
    File? pango_lineage_report
    String organism = "sars-cov-2"
    Int disk_size = 100
  }
  command <<<
    # when running on terra, comment out all input_table mentions
    #python3 /scripts/export_large_tsv/export_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --entity_type ~{table_name} --tsv_filename ~{table_name}-data.tsv
    
    # when running locally, use the input_table in place of downloading from Terra
    cp ~{input_table} ~{table_name}-data.tsv

    echo "DEBUG: Now entering Python block to perform parsing of data for populating sample-level table"

    python3 <<CODE 
    import pandas as pd 
    import numpy as np 
    import csv
    import re
    import os 
    import sys

    # read exported Terra table into pandas
    tablename = "~{table_name}-data.tsv" 
    table = pd.read_csv(tablename, delimiter='\t', header=0, dtype={"~{table_name}_id": 'str'}) # ensure sample_id is always a string

    # extract the samples for upload from the entire table
    table = table[table["~{table_name}_id"].isin("~{sep='*' samplenames}".split("*"))]

    # set all column headers to lowercase 
    table.columns = table.columns.str.lower()

    # parse nextclade into dataframe if path exists
    if os.path.exists("~{nextclade_tsv}"):
      print("Nextclade tsv exists, parsing into dataframe")
      nextclade = pd.read_csv("~{nextclade_tsv}", delimiter='\t')
    
    # parse pango lineage report into dataframe if path exists
    if os.path.exists("~{pango_lineage_report}"):
      print("Pango lineage report exists, parsing into dataframe")
      pango_lineage_report = pd.read_csv("~{pango_lineage_report}", delimiter=',')
    
    # set required and optional metadata fields based on the organism type
    if ("~{organism}" == "sars-cov-2"):
      # TODO: Assuming nextclade and pango lineage reports exist for all SARS-CoV-2 samples
      print("Organism is SARS-CoV-2, parsing nexclade and pango lineage reports") # missing VADR
      
      for index, row in table.iterrows():
        # TODO: assuming assembly fastas are in assembly_fasta row - should change this to be more flexible
        assembly_name = row["assembly_fasta"].split('/')[-1].split('.')[0]
        
        # NEXTCLADE PARSING BLOCK
        # TODO - missing output files: nextclade_json, auspice_json, nextclade_tsv
        #   One file for all samples, would have to split into multiple files, upload to bucket
        # TODO - missing output strings: nextclade_version, nextclade_docker, nextclade_ds_tag

        # parse nextclade_clade from nextclade tsv
        nc_clade = nextclade.loc[nextclade["seqName"] == assembly_name]["clade_nextstrain"].item()
        who_clade = nextclade.loc[nextclade["seqName"] == assembly_name]["clade_who"].item()
        if (nc_clade != who_clade) and (nc_clade != '') and (who_clade != ''):
          nc_clade = nc_clade + " (" + who_clade + ")"
        if nc_clade == '':
          nc_clade = 'NA'
        # replace nextclade value in datatable if exists, if not, create it
        if "nextclade_clade" not in table.columns:
          table["nextclade_clade"] = ""
        table["nextclade_clade"][index] = nc_clade

        # parse nexclade nextclade_aa_subs from nextclade tsv
        nc_aa_subs = nextclade.loc[nextclade["seqName"] == assembly_name]["aaSubstitutions"].item()
        if nc_aa_subs == '':
          nc_aa_subs = 'NA'
        # replace nextclade value in datatable if exists, if not, create it
        if "nextclade_aa_subs" not in table.columns:
          table["nextclade_aa_subs"] = ""
        table["nextclade_aa_subs"][index] = nc_aa_subs

        # parse nexclade nextclade_aa_dels from nextclade tsv
        nc_aa_dels = nextclade.loc[nextclade["seqName"] == assembly_name]["aaDeletions"].item()
        if nc_aa_dels == '':
          nc_aa_dels = 'NA'
        # replace nextclade value in datatable if exists, if not, create it
        if "nextclade_aa_dels" not in table.columns:
          table["nextclade_aa_dels"] = ""
        table["nextclade_aa_dels"][index] = nc_aa_dels

        # parse nexclade nextclade_lineage from nextclade tsv
        try:
          nc_lineage = nextclade.loc[nextclade["seqName"] == assembly_name]["lineage"].item()
        except KeyError:
          nc_lineage = ''
        if nc_lineage == '':
          nc_lineage = 'NA'
        # replace nextclade value in datatable if exists, if not, create it
        if "nextclade_lineage" not in table.columns:
          table["nextclade_lineage"] = ""
        table["nextclade_lineage"][index] = nc_lineage

        # PANGO LINEAGE PARSING BLOCK
        # TODO: missing output files: pango_lineage_report
        # TODO: missing output strings: pangolin_docker

        # parse pangolin_version from pango lineage report
        pangolin_version = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["pangolin_version"].item()
        version = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["version"].item()
        assignment_version = "pangolin {pangolin_version}; {version}".format(pangolin_version=pangolin_version, version=version)
        # replace pangolin_version in datatable if exists, if not, create it
        if "pangolin_version" not in table.columns:
          table["pangolin_version"] = ""
        table["pangolin_version"][index] = assignment_version

        # parse pango_lineage from pango lineage report
        pango_lineage = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["lineage"].item()
        if "pango_lineage" not in table.columns:
          table["pango_lineage"] = ""
        table["pango_lineage"][index] = pango_lineage

        # parse pango_lineage_expanded from pango lineage report
        try:
          pango_lineage_expanded = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["expanded_lineage"].item()
        except KeyError:
            pango_lineage_expanded = ""
        if "pango_lineage_expanded" not in table.columns:
          table["pango_lineage_expanded"] = ""
        table["pango_lineage_expanded"][index] = pango_lineage_expanded

        # parse pangolin_conflicts from pango lineage report
        pangolin_conflicts = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["conflict"].item()
        if "pangolin_conflicts" not in table.columns:
          table["pangolin_conflicts"] = ""
        table["pangolin_conflicts"][index] = pangolin_conflicts

        # parse pangolin_notes from pango lineage report
        pangolin_notes = pango_lineage_report.loc[pango_lineage_report["taxon"] == assembly_name]["note"].item()
        if "pangolin_notes" not in table.columns:
          table["pangolin_notes"] = ""
        table["pangolin_notes"][index] = pangolin_notes
    
    # TODO: drop assembly_fasta column? What about other columns that are not required for upload?
    table.to_csv("TERRA_TABLE_TEMP.tsv", sep='\t', index=False)

    CODE

    # upload results to terra databable 
    #python3 scripts/import_large_tsv/import_large_tsv.py --project "~{project_name}" --workspace "~{workspace_name}" --tsv <path_to_tsv_to_upload> 

    echo "DEBUG: upload to terra table complete"
  >>>
  output {
    File terra_table = "TERRA_TABLE_TEMP.tsv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}