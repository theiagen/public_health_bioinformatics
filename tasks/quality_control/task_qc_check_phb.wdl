version 1.0

task qc_check_phb {
  input {
    File? qc_check_table
    String? gambit_predicted_taxon
    String? expected_taxon
    Int? assembly_length
    Float? est_coverage_raw
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

        if ("assembly_length_min" in qc_check_metrics) and ("assembly_length_max" in qc_check_metrics):
          if (~{assembly_length}):
            qc_note, qc_status = compare(qc_note, "assembly_length", int(~{assembly_length}), ">=", int(taxon_df["assembly_length_min"][0]), "<=", int(taxon_df["assembly_length_max"][0]))
            qc_check_metrics.remove("assembly_length_min")
            qc_check_metrics.remove("assembly_length_max")

        if ("est_coverage_raw" in qc_check_metrics): # if this var is in the qc_check_metrics,
          if (~{est_coverage_raw}): # if est_coverage_raw variable exists,
            qc_note, qc_status = compare(qc_note, "est_coverage_raw", float(~{est_coverage_raw}), ">=", float(taxon_df["est_coverage_raw"][0]))
            qc_check_metrics.remove("est_coverage_raw")

        # add more measures here

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