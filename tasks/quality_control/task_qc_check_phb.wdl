version 1.0

task qc_check_phb {
  input {
    # workflow agnotic inputs
    File? qc_check_table
    String? gambit_predicted_taxon
    String? expected_taxon
    # theiaprok and theiaeuk inputs
    Float? r1_mean_q_raw
    Float? r2_mean_q_raw
    Float? combined_mean_q_raw
    Float? r1_mean_readlength_raw
    Float? r2_mean_readlength_raw    
    Float? combined_mean_readlength_raw 
    Float? r1_mean_q_clean
    Float? r2_mean_q_clean
    Float? combined_mean_q_clean
    Float? r1_mean_readlength_clean
    Float? r2_mean_readlength_clean    
    Float? combined_mean_readlength_clean 
    Float? est_coverage_raw
    Float? est_coverage_clean 
    Int? assembly_length
    Int? number_contigs 
    Int? n50_value 
    String? busco_results
    # theiaprok only inputs
    String? midas_secondary_genus_abundance 
    Float? ani_highest_percent 
    Float? ani_highest_percent_bases_aligned
    # theiacov inputs
    Int? num_reads_raw1
    Int? num_reads_raw2
    Int? num_reads_clean1
    Int? num_reads_clean2
    Float? kraken_human 
    Float? kraken_human_dehosted
    Float? kraken_sc2
    Float? kraken_sc2_dehosted
    Float? kraken_target_org
    Float? kraken_target_org_dehosted
    String? meanbaseq_trim
    Float? assembly_mean_coverage
    Int? number_N
    Int? number_Degenerate
    Int? assembly_length_unambiguous
    Float? percent_reference_coverage
    Float? sc2_s_gene_mean_coverage
    Float? sc2_s_gene_percent_coverage
    String? vadr_num_alerts
    Int disk_size = 100
  }
  command <<<
    python3 <<CODE
    import csv
    import pandas as pd
    import numpy as np

    # set a function to compare the input to a standard value
    # qc_note: the notes regarding the qc_check to be appended to
    # input: input value to examine (already cast to intended type)
    # expectation: should this input be >, >=, =, <, <= to the standard
    # standard: the value to compare the input to
    # second_expectation: should this input be >, >=, =, <, <= to the second_standard
    # second_standard: in case there are two values to compare (e.g., upper/lower boundary)
    def compare(qc_note, variable_name, input_value, expectation, standard, second_expectation=None, second_standard=None):
      # create empty variable to return
      qc_status = ""

      # check if input value exists
      if (input_value): 
        # perform check on every possible operator
        if expectation == ">":
          if (input_value > standard):
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was greater than the threshold of " + str(standard))
            if (second_standard): # if a second boundary exists, recursively add to the qc_note line
              qc_note += compare(qc_note, variable_name, input_value, second_expectation, second_standard)[0]
          else:
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was not greater than the threshold of " + str(standard))
            qc_note += variable_name + " (" + str(input_value) + ") was less than or equal to the minimum threshold of " + str(standard) + "; "

        elif expectation == ">=":
          if (input_value >= standard):
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was greater than or equal to the threshold of " + str(standard))
            if (second_standard): # if a second boundary exists, recursively add to the qc_note line
              qc_note += compare(qc_note, variable_name, input_value, second_expectation, second_standard)[0]
          else:
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was not greater than or equal to the threshold of " + str(standard))
            qc_note += variable_name + " (" + str(input_value) + ") was less than the minimum threshold of " + str(standard) + "; "

        elif expectation == "=":
          if (input_value == standard):
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was equal to the threshold of " + str(standard))
            if (second_standard): # if a second boundary exists, recursively add to the qc_note line
              qc_note += compare(qc_note, variable_name, input_value, second_expectation, second_standard)[0]
          else:
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was not equal to the threshold of " + str(standard))
            qc_note += variable_name + " (" + str(input_value) + ") was not equal to the threshold of " + str(standard) + "; "
            
        elif expectation == "<":
          if (input_value < standard):
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was less than the threshold of " + str(standard))
            if (second_standard): # if a second boundary exists, recursively add to the qc_note line
              qc_note += compare(qc_note, variable_name, input_value, second_expectation, second_standard)[0]
          else:
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was not less than the threshold of " + str(standard))
            qc_note += variable_name + " (" + str(input_value) + ") was greater than or equal to the maximum threshold of " + str(standard) + "; "
            
        elif expectation == "<=":
          if (input_value <= standard):
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was less than or equal to the threshold of " + str(standard))
            if (second_standard): # if a second boundary exists, recursively add to the qc_note line
              qc_note += compare(qc_note, variable_name, input_value, second_expectation, second_standard)[0]
          else:
            print("DEBUG: " + variable_name + " (" + str(input_value) + ") was not less or equal to the threshold of " + str(standard))
            qc_note += variable_name + " (" + str(input_value) + ") was greater than the maximum threshold of " + str(standard) + "; "

      # if the qc_note has a value, then it has failed a check
      if (len(qc_note) > 0):
        qc_status = "QC_ALERT"

      return qc_note, qc_status

    # create two empty variables for results to be added to
    qc_status = ""
    qc_note = ""
    
    # import the qc_check_table into a pandas data frame
    qc_check_df = pd.read_csv("~{qc_check_table}", sep = '\t', index_col = "taxon")
    
    # extract the list of taxon to examine
    qc_check_taxa = qc_check_df.index.values.tolist()

    # preferentially use a user-provided expected taxon; otherwise, use the gambit predicted taxon
    expected_taxon = "~{expected_taxon}"
    gambit_predicted_taxon = "~{gambit_predicted_taxon}"

    # check if an expected_taxon or gambit_predicted taxon was detected
    if (expected_taxon):
      qc_taxon = expected_taxon.replace(" ", "_")
      print("DEBUG: User-provided expected_taxon was found and will be used for QC check")
    elif (gambit_predicted_taxon):
      qc_taxon = gambit_predicted_taxon.replace(" ", "_")
      print("DEBUG: No user-provided expected_taxon found; will be using gambit_predicted_taxon for QC check")
    else:
      qc_status = "QC_NA"
      qc_note = "No expected_taxon or gambit_predicted_taxon was found; qc_check could not proceed"
    
    # if a qc_taxon was generated, check to see if that taxon is in the qc_check_table
    if (qc_status != "QC_NA"):
      matching_taxon = []
      for taxa in qc_check_taxa:
        if taxa in qc_taxon:
          matching_taxon.append(taxa)
          print("DEBUG: the qc_taxon was found, adding to matching_taxon list")
      if len(matching_taxon) < 1: # if no taxon were found
        qc_status = "QC_NA"
        qc_note = "No matching taxon were found in the qc_check table; qc_check could not proceed"
      elif len(matching_taxon) > 1: # if more than one taxon were found
        qc_status = "QC_NA"
        qc_note = "Multiple matching taxa were detected in the qc_check_table; qc_check could not proceed"
      else: # only one taxon was found
        taxon_df = qc_check_df.loc[[matching_taxon[0]]] # extract only the values for that taxa
        print(f"DEBUG: exactly one matching taxon was found: {matching_taxon}, proceeding with qc_check")

        # remove columns where all values are null
        taxon_df = taxon_df.replace(r'^\s*$', np.nan, regex=True)
        taxon_df = taxon_df.dropna(how='all', axis=1)
        #print(taxon_df)
        
      # perform qc_check on any metrics in the qc_check_table
      if (qc_status != "QC_NA"):
        qc_check_metrics = taxon_df.columns.values.tolist()
        #print(qc_check_metrics)

        if ("r1_mean_q_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r1_mean_q_raw}): # if r1_mean_q_raw variable exists,
            qc_note, qc_status = compare(qc_note, "r1_mean_q_raw", float(~{r1_mean_q_raw}), ">=", float(taxon_df["r1_mean_q_raw"][0]))
            qc_check_metrics.remove("r1_mean_q_raw")

        if ("r2_mean_q_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r2_mean_q_raw}): # if r2_mean_q_raw variable exists,
            qc_note, qc_status = compare(qc_note, "r2_mean_q_raw", float(~{r2_mean_q_raw}), ">=", float(taxon_df["r2_mean_q_raw"][0]))
            qc_check_metrics.remove("r2_mean_q_raw")

        if ("combined_mean_q_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{combined_mean_q_raw}): # if combined_mean_q_raw variable exists,
            qc_note, qc_status = compare(qc_note, "combined_mean_q_raw", float(~{combined_mean_q_raw}), ">=", float(taxon_df["combined_mean_q_raw"][0]))
            qc_check_metrics.remove("combined_mean_q_raw")

        if ("r1_mean_readlength_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r1_mean_readlength_raw}): # if r1_mean_readlength_raw variable exists,
            qc_note, qc_status = compare(qc_note, "r1_mean_readlength_raw", float(~{r1_mean_readlength_raw}), ">=", float(taxon_df["r1_mean_readlength_raw"][0]))
            qc_check_metrics.remove("r1_mean_readlength_raw")

        if ("r2_mean_readlength_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r2_mean_readlength_raw}): # if r2_mean_readlength_raw variable exists,
            qc_note, qc_status = compare(qc_note, "r2_mean_readlength_raw", float(~{r2_mean_readlength_raw}), ">=", float(taxon_df["r2_mean_readlength_raw"][0]))
            qc_check_metrics.remove("r2_mean_readlength_raw")

        if ("combined_mean_readlength_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{combined_mean_readlength_raw}): # if combined_mean_readlength_raw variable exists,
            qc_note, qc_status = compare(qc_note, "combined_mean_readlength_raw", float(~{combined_mean_readlength_raw}), ">=", float(taxon_df["combined_mean_readlength_raw"][0]))
            qc_check_metrics.remove("combined_mean_readlength_raw")

        if ("r1_mean_q_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r1_mean_q_clean}): # if r1_mean_q_clean variable exists,
            qc_note, qc_status = compare(qc_note, "r1_mean_q_clean", float(~{r1_mean_q_clean}), ">=", float(taxon_df["r1_mean_q_clean"][0]))
            qc_check_metrics.remove("r1_mean_q_clean")

        if ("r2_mean_q_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r2_mean_q_clean}): # if r2_mean_q_clean variable exists,
            qc_note, qc_status = compare(qc_note, "r2_mean_q_clean", float(~{r2_mean_q_clean}), ">=", float(taxon_df["r2_mean_q_clean"][0]))
            qc_check_metrics.remove("r2_mean_q_clean")

        if ("combined_mean_q_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{combined_mean_q_clean}): # if combined_mean_q_clean variable exists,
            qc_note, qc_status = compare(qc_note, "combined_mean_q_clean", float(~{combined_mean_q_clean}), ">=", float(taxon_df["combined_mean_q_clean"][0]))
            qc_check_metrics.remove("combined_mean_q_clean")

        if ("r1_mean_readlength_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r1_mean_readlength_clean}): # if r1_mean_readlength_clean variable exists,
            qc_note, qc_status = compare(qc_note, "r1_mean_readlength_clean", float(~{r1_mean_readlength_clean}), ">=", float(taxon_df["r1_mean_readlength_clean"][0]))
            qc_check_metrics.remove("r1_mean_readlength_clean")

        if ("r2_mean_readlength_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{r2_mean_readlength_clean}): # if r2_mean_readlength_clean variable exists,
            qc_note, qc_status = compare(qc_note, "r2_mean_readlength_clean", float(~{r2_mean_readlength_clean}), ">=", float(taxon_df["r2_mean_readlength_clean"][0]))
            qc_check_metrics.remove("r2_mean_readlength_clean")

        if ("combined_mean_readlength_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{combined_mean_readlength_clean}): # if combined_mean_readlength_clean variable exists,
            qc_note, qc_status = compare(qc_note, "combined_mean_readlength_clean", float(~{combined_mean_readlength_clean}), ">=", float(taxon_df["combined_mean_readlength_clean"][0]))
            qc_check_metrics.remove("combined_mean_readlength_clean")

        if ("est_coverage_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{est_coverage_raw}): # if est_coverage_raw variable exists,
            qc_note, qc_status = compare(qc_note, "est_coverage_raw", float(~{est_coverage_raw}), ">=", float(taxon_df["est_coverage_raw"][0]))
            qc_check_metrics.remove("est_coverage_raw")

        if ("est_coverage_clean" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{est_coverage_clean}): # if est_coverage_clean variable exists,
            qc_note, qc_status = compare(qc_note, "est_coverage_clean", float(~{est_coverage_clean}), ">=", float(taxon_df["est_coverage_clean"][0]))
            qc_check_metrics.remove("est_coverage_clean")

        if ("midas_secondary_genus_abundance" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{midas_secondary_genus_abundance}): # if midas_secondary_genus_abundance variable exists,
            qc_note, qc_status = compare(qc_note, "midas_secondary_genus_abundance", float(~{midas_secondary_genus_abundance}), "<", float(taxon_df["midas_secondary_genus_abundance"][0]))
            qc_check_metrics.remove("midas_secondary_genus_abundance")

        if ("assembly_length_min" in qc_check_metrics) and ("assembly_length_max" in qc_check_metrics):
          if (~{assembly_length}):
            qc_note, qc_status = compare(qc_note, "assembly_length", int(~{assembly_length}), ">=", int(taxon_df["assembly_length_min"][0]), "<=", int(taxon_df["assembly_length_max"][0]))
            qc_check_metrics.remove("assembly_length_min")
            qc_check_metrics.remove("assembly_length_max")     

        if ("number_contigs" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{number_contigs}): # if number_contigs variable exists,
            qc_note, qc_status = compare(qc_note, "number_contigs", int(~{number_contigs}), "<=", int(taxon_df["number_contigs"][0]))
            qc_check_metrics.remove("number_contigs")

        if ("n50_value" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{n50_value}): # if n50_value variable exists,
            qc_note, qc_status = compare(qc_note, "n50_value", int(~{n50_value}), ">=", int(taxon_df["n50_value"][0]))
            qc_check_metrics.remove("n50_value")

        if ("ani_highest_percent" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{ani_highest_percent}): # if ani_highest_percent variable exists,
            qc_note, qc_status = compare(qc_note, "ani_highest_percent", float(~{ani_highest_percent}), ">=", float(taxon_df["ani_highest_percent"][0]))
            qc_check_metrics.remove("ani_highest_percent")

        if ("ani_highest_percent_bases_aligned" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{ani_highest_percent_bases_aligned}): # if ani_highest_percent_bases_aligned variable exists,
            qc_note, qc_status = compare(qc_note, "ani_highest_percent_bases_aligned", float(~{ani_highest_percent_bases_aligned}), ">=", float(taxon_df["ani_highest_percent_bases_aligned"][0]))
            qc_check_metrics.remove("ani_highest_percent_bases_aligned")

        if ("busco_completeness" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if ("~{busco_results}"): # if busco_results variable exists, quotes are necessary because of invalid syntax error
            # parse busco_results string into busco completeness value only
            busco_completeness = "~{busco_results}".split("C:")[1].split("%")[0]
            qc_note, qc_status = compare(qc_note, "busco_completeness", float(busco_completeness), ">=", float(taxon_df["busco_completeness"][0]))
            qc_check_metrics.remove("busco_completeness")
        
        if ("num_reads_raw1" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{num_reads_raw1}): # if num_reads_raw1 variable exists,
            qc_note, qc_status = compare(qc_note, "num_reads_raw1", float(~{num_reads_raw1}), ">=", float(taxon_df["num_reads_raw1"][0]))
            qc_check_metrics.remove("num_reads_raw1")

        if ("num_reads_raw2" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{num_reads_raw2}): # if num_reads_raw2 variable exists,
            qc_note, qc_status = compare(qc_note, "num_reads_raw2", int(~{num_reads_raw2}), ">=", int(taxon_df["num_reads_raw2"][0]))
            qc_check_metrics.remove("num_reads_raw2") 

        if ("num_reads_clean1" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{num_reads_clean1}): # if num_reads_clean1 variable exists,
            qc_note, qc_status = compare(qc_note, "num_reads_clean1", int(~{num_reads_clean1}), ">=", int(taxon_df["num_reads_clean1"][0]))
            qc_check_metrics.remove("num_reads_clean1")

        if ("num_reads_clean2" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{num_reads_clean2}): # if num_reads_clean2 variable exists,
            qc_note, qc_status = compare(qc_note, "num_reads_clean2", int(~{num_reads_clean2}), ">=", int(taxon_df["num_reads_clean2"][0]))
            qc_check_metrics.remove("num_reads_clean2")

        if ("kraken_human" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_human}): # if kraken_human variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_human", float(~{kraken_human}), "<=", float(taxon_df["kraken_human"][0]))
            qc_check_metrics.remove("kraken_human")

        if ("kraken_human_dehosted" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_human_dehosted}): # if kraken_human_dehosted variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_human_dehosted", float(~{kraken_human_dehosted}), "<=", float(taxon_df["kraken_human_dehosted"][0]))
            qc_check_metrics.remove("kraken_human_dehosted")

        if ("kraken_sc2" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_sc2}): # if kraken_sc2 variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_sc2", float(~{kraken_sc2}), ">=", float(taxon_df["kraken_sc2"][0]))
            qc_check_metrics.remove("kraken_sc2") 

        if ("kraken_sc2_dehosted" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_sc2_dehosted}): # if kraken_sc2_dehosted variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_sc2_dehosted", float(~{kraken_sc2_dehosted}), ">=", float(taxon_df["kraken_sc2_dehosted"][0]))
            qc_check_metrics.remove("kraken_sc2_dehosted")   

        if ("kraken_target_org" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_target_org}): # if kraken_target_org variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_target_org", float(~{kraken_target_org}), ">=", float(taxon_df["kraken_target_org"][0]))
            qc_check_metrics.remove("kraken_target_org") 

        if ("kraken_target_org_dehosted" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{kraken_target_org_dehosted}): # if kraken_target_org_dehosted variable exists,
            qc_note, qc_status = compare(qc_note, "kraken_target_org_dehosted", float(~{kraken_target_org_dehosted}), ">=", float(taxon_df["kraken_target_org_dehosted"][0]))
            qc_check_metrics.remove("kraken_target_org_dehosted")   

        if ("meanbaseq_trim" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{meanbaseq_trim}): # if meanbaseq_trim variable exists,
            qc_note, qc_status = compare(qc_note, "meanbaseq_trim", float(~{meanbaseq_trim}), ">=", float(taxon_df["meanbaseq_trim"][0]))
            qc_check_metrics.remove("meanbaseq_trim")  

        if ("assembly_mean_coverage" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{assembly_mean_coverage}): # if assembly_mean_coverage variable exists,
            qc_note, qc_status = compare(qc_note, "assembly_mean_coverage", float(~{assembly_mean_coverage}), ">=", float(taxon_df["assembly_mean_coverage"][0]))
            qc_check_metrics.remove("assembly_mean_coverage")  

        if ("number_N" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{number_N}): # if number_N variable exists,
            qc_note, qc_status = compare(qc_note, "number_N", int(~{number_N}), "<=", int(taxon_df["number_N"][0]))
            qc_check_metrics.remove("number_N") 

        if ("number_Degenerate" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{number_Degenerate}): # if number_Degenerate variable exists,
            qc_note, qc_status = compare(qc_note, "number_Degenerate", int(~{number_Degenerate}), "<=", int(taxon_df["number_Degenerate"][0]))
            qc_check_metrics.remove("number_Degenerate") 

        if ("assembly_length_unambiguous_min" in qc_check_metrics) and ("assembly_length_unambiguous_max" in qc_check_metrics):
          if (~{assembly_length}):
            qc_note, qc_status = compare(qc_note, "assembly_length_unambiguous", int(~{assembly_length_unambiguous}), ">=", int(taxon_df["assembly_length_unambiguous_min"][0]), "<=", int(taxon_df["assembly_length_unambiguous_max"][0]))
            qc_check_metrics.remove("assembly_length_unambiguous_min")
            qc_check_metrics.remove("assembly_length_unambiguous_max")     

        if ("percent_reference_coverage" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{percent_reference_coverage}): # if percent_reference_coverage variable exists,
            qc_note, qc_status = compare(qc_note, "percent_reference_coverage", float(~{percent_reference_coverage}), ">=", float(taxon_df["percent_reference_coverage"][0]))
            qc_check_metrics.remove("percent_reference_coverage") 

        if ("sc2_s_gene_mean_coverage" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{sc2_s_gene_mean_coverage}): # if sc2_s_gene_mean_coverage variable exists,
            qc_note, qc_status = compare(qc_note, "sc2_s_gene_mean_coverage", float(~{sc2_s_gene_mean_coverage}), ">=", float(taxon_df["sc2_s_gene_mean_coverage"][0]))
            qc_check_metrics.remove("sc2_s_gene_mean_coverage") 

        if ("sc2_s_gene_percent_coverage" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{sc2_s_gene_percent_coverage}): # if sc2_s_gene_percent_coverage variable exists,
            qc_note, qc_status = compare(qc_note, "sc2_s_gene_percent_coverage", float(~{sc2_s_gene_percent_coverage}), ">=", float(taxon_df["sc2_s_gene_percent_coverage"][0]))
            qc_check_metrics.remove("sc2_s_gene_percent_coverage") 

        if ("vadr_num_alerts" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{vadr_num_alerts}): # if vadr_num_alerts variable exists,
            qc_note, qc_status = compare(qc_note, "vadr_num_alerts", int(~{vadr_num_alerts}), ">=", int(taxon_df["vadr_num_alerts"][0]))
            qc_check_metrics.remove("vadr_num_alerts") 

        if (len(qc_check_metrics) > 0):
          qc_status = "QC_ALERT"
          qc_note += f"one or more columns in qc_check_metrics was missing a required input: {qc_check_metrics}"

      if (qc_status != "QC_NA") and (qc_status != "QC_ALERT"):
        qc_check = "QC_PASS"
      else:
        qc_check = qc_status + ": " + qc_note
        qc_check = qc_check.rstrip('; ')

      with open("QC_CHECK", 'wt') as out:
        out.write(qc_check)
      
    CODE

  >>>
  output {
    String qc_check = read_string("QC_CHECK")
    File? qc_standard = qc_check_table
  }
  runtime {
    docker: "quay.io/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
   # maxRetries: 3
    preemptible: 0
  }
}