version 1.0

task qc_check_phb {
  input {
    # core inputs
    File? qc_check_table
    # {qc_metric: [value, type, operator, use_exception]}
    Map[String, String?] qc_check_inputs
    File? irma_qc_table

    String? expected_taxon
    String? gambit_predicted_taxon

    Int disk_size = 100
    Int memory = 8
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2024-08-27"
  }
  File qc_check_input_json = write_json(qc_check_inputs)
  command <<<
    python3 <<CODE
    import csv
    import json
    import operator
    import pandas as pd
    import numpy as np

    # define the qc_check_criteria metadata (lower-cased keys for case-insensitive parsing)
    # errors SHOULD occur if key is not present here to flag devs to add new checks
    # "_min" and "_max" suffixes are allowed in qc_check_inputs to override operator
    qc_check_criteria = {
      "ani_highest_percent": {"type": float, "operator": operator.ge},
      "ani_highest_percent_bases_aligned": {"type": float, "operator": operator.ge},
      "assembly_length": {"type": int, "operator": operator.ge},
      "assembly_length_unambiguous": {"type": int, "operator": operator.ge},
      "assembly_mean_coverage": {"type": float, "operator": operator.ge},
      "busco_completeness": {"type": float, "operator": operator.ge},
      "combined_mean_q_clean": {"type": float, "operator": operator.ge},
      "combined_mean_q_raw": {"type": float, "operator": operator.ge},
      "combined_mean_readlength_clean": {"type": float, "operator": operator.ge},
      "combined_mean_readlength_raw": {"type": float, "operator": operator.ge},
      "est_coverage_clean": {"type": float, "operator": operator.ge},
      "est_coverage_raw": {"type": float, "operator": operator.ge},
      "kraken_human": {"type": float, "operator": operator.le},
      "kraken_human_dehosted": {"type": float, "operator": operator.le},
      "meanbaseq_trim": {"type": float, "operator": operator.ge},
      "midas_secondary_genus_abundance": {"type": float, "operator": operator.lt},
      "midas_secondary_genus_coverage": {"type": float, "operator": operator.lt},
      "n50_value": {"type": int, "operator": operator.ge},
      "num_reads_cleaned1": {"type": int, "operator": operator.ge},
      "num_reads_cleaned2": {"type": int, "operator": operator.ge},
      "num_reads_raw1": {"type": int, "operator": operator.ge},
      "num_reads_raw2": {"type": int, "operator": operator.ge},
      "number_contigs": {"type": int, "operator": operator.le},
      "number_degenerate": {"type": int, "operator": operator.le},
      "number_n": {"type": int, "operator": operator.le},
      "percent_reference_coverage": {"type": float, "operator": operator.ge},
      "quast_gc_percent": {"type": float, "operator": operator.ge},
      "r1_mean_q_clean": {"type": float, "operator": operator.ge},
      "r1_mean_q_raw": {"type": float, "operator": operator.ge},
      "r1_mean_readlength_clean": {"type": float, "operator": operator.ge},
      "r1_mean_readlength_raw": {"type": float, "operator": operator.ge},
      "r2_mean_q_clean": {"type": float, "operator": operator.ge},
      "r2_mean_q_raw": {"type": float, "operator": operator.ge},
      "r2_mean_readlength_clean": {"type": float, "operator": operator.ge},
      "r2_mean_readlength_raw": {"type": float, "operator": operator.ge},
      "sc2_s_gene_mean_coverage": {"type": float, "operator": operator.ge},
      "sc2_s_gene_percent_coverage": {"type": float, "operator": operator.ge},
      "vadr_num_alerts": {"type": int, "operator": operator.le}
    }
    # create a map for segment-based checks
    segment_check_criteria = {
      "median_coverage": {"type": float, "operator": operator.ge, "col": "Median Coverage"},
      "num_minor_snv": {"type": int, "operator": operator.le, "col": "Count of Minor SNVs (AF >= 0.05)"},
      "percent_reference_coverage": {"type": float, "operator": operator.ge, "col": "% Reference Covered"}
    }

    # set a function to compare the input to a standard value
    # qc_note: the notes regarding the qc_check to be appended to
    # input: input value to examine (already cast to intended type)
    # expectation: should this input be >, >=, =, <, <= to the standard
    # standard: the value to compare the input to
    def compare(qc_note, variable_name, input_value, operator, standard):
      # create empty variable to return
      qc_status = ""

      # operator str reporting map
      op2str = {
        operator.gt: "greater than",
        operator.ge: "greater than or equal to",
        operator.eq: "equal to",
        operator.lt: "less than",
        operator.le: "less than or equal to"
      }
      # operator inverse str reporting map
      op2inv_str = {
        operator.gt: "less than or equal to the minimum",
        operator.ge: "less than the minimum",
        operator.eq: "not equal to",
        operator.lt: "greater than or equal to the maximum",
        operator.le: "greater than the maximum"
      }

      # check if input value exists
      if input_value is not None and not pd.isnull(input_value): 
        # perform check on every possible operator
        if operator(input_value, standard):
          print(f"DEBUG: {variable_name} ({input_value}) was {op2str[operator]} the threshold of {standard}")
        else:
          print(f"DEBUG: {variable_name} ({input_value}) was not {op2str[operator]} the threshold of {standard}")
          qc_note += f"{variable_name} ({input_value}) was {op2inv_str[operator]} threshold of {standard}; "

      # if the qc_note has a value, then it has failed a check
      if qc_note:
        qc_status = "QC_ALERT"

      return qc_note, qc_status

    # create two empty variables for results to be added to
    qc_status = ""
    qc_note = ""
    
    # import the qc_check_table into a pandas data frame
    qc_check_df = pd.read_csv("~{qc_check_table}", sep = '\t', index_col = "taxon")

    # lower-case for case-insensitive parsing
    qc_check_df.columns = qc_check_df.columns.str.lower()
    
    # import and clean the qc_check_criteria json generated by WDL
    with open("~{qc_check_input_json}", "r") as f:
      qc_check_inputs_dirty = json.load(f)
    qc_check_inputs = {}
    for wdl_dict in qc_check_inputs_dirty:
      key = wdl_dict["left"]
      obs_value = wdl_dict["right"]
      # skip NoneType (WDL null) values
      if obs_value is not None:
        qc_check_inputs[key.lower()] = obs_value

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
        
      # perform qc_check on any metrics in the qc_check_table
      if (qc_status != "QC_NA"):
        qc_check_metrics = taxon_df.columns.values.tolist()
        print(f"DEBUG: Found qc_check_metrics: {qc_check_metrics}")

        # iterate through standard checks first
        for metric in sorted(qc_check_metrics):
          # allow "min" and "max" suffixes to be used in qc_check_criteria
          if metric.endswith("_min"):
            obs_operator = operator.ge
            base_metric = metric[:-4]
          elif metric.endswith("_max"):
            obs_operator = operator.le
            base_metric = metric[:-4]
          else:
            obs_operator = None
            base_metric = metric
          if base_metric in qc_check_inputs:
            obs_val = qc_check_inputs[base_metric]
            if base_metric not in qc_check_criteria:
              raise ValueError(f"qc_check_criteria does not contain an entry for {base_metric}; please add it to the qc_check_criteria")
            else:
              val_type = qc_check_criteria[base_metric]["type"]
              # yield to the operator from "min" / "max" suffixes if present
              if not obs_operator:
                obs_operator = qc_check_criteria[base_metric]["operator"]
            # special exemption we want to avoid as much as possible
            if base_metric == "busco_completeness":
              busco_completeness = qc_check_criteria[base_metric][0].split("C:")[1].split("%")[0]
              qc_note, qc_status = compare(qc_note, metric, float(busco_completeness), obs_operator, float(taxon_df[metric][0]))
            else:
              # catch casting errors, particularly for VADR
              try:
                qc_note, qc_status = compare(qc_note, metric, val_type(obs_val), obs_operator, val_type(taxon_df[metric][0]))
              except:
                qc_note += f"{base_metric} ({obs_val}) could not be cast to {val_type_str}; "
                qc_status = "QC_ALERT"
            qc_check_metrics.remove(metric)

        # check segment-level qc metrics via the irma_qc_table if it exists
        if "~{irma_qc_table}":
          irma_qc_table = pd.read_csv("~{irma_qc_table}", sep = '\t', index_col = "Sample")
          to_rm = set()
          # report segment failures on a metric-by-metric basis, but we iterate segment-by-segment
          seg2qc_note = {}
          # iterate through segments
          for index, row in irma_qc_table.iterrows():
            # extract the segment name
            if '_' in row['Reference']:
              segment_name = row['Reference'].lower()[:row['Reference'].rfind('_')]
            else:
              segment_name = row['Reference'].lower()
            print("DEBUG: Checking QC metrics for segment: " + segment_name)
            seg2qc_note[segment_name] = {}
            # iterate through the segment check criteria
            for var, metadata in segment_check_criteria.items():
              val_type = metadata["type"]
              obs_operator = metadata["operator"]
              col_name = metadata["col"]
              obs_val = row[col_name]
              full_var_name = f"{segment_name}_{var}"
              gen_var_name = f"segment_{var}"
              seg_qc_note = ""
              # prioritize segment-specific thresholds
              if full_var_name in qc_check_metrics:
                try:
                  seg_qc_note, seg_qc_status = compare(seg_qc_note, full_var_name, val_type(obs_val), obs_operator, val_type(taxon_df[full_var_name][0]))
                # cast to a float if initial cast fails, usually due to NaNs
                except ValueError:
                  seg_qc_note, seg_qc_status = compare(seg_qc_note, full_var_name, float(obs_val), obs_operator, val_type(taxon_df[full_var_name][0]))
                to_rm.add(full_var_name)
              elif gen_var_name in qc_check_metrics:
                try:
                  seg_qc_note, seg_qc_status = compare(seg_qc_note, gen_var_name, val_type(obs_val), obs_operator, val_type(taxon_df[gen_var_name][0]))
                except ValueError:
                  seg_qc_note, seg_qc_status = compare(seg_qc_note, gen_var_name, float(obs_val), obs_operator, val_type(taxon_df[gen_var_name][0]))
                to_rm.add(gen_var_name)
              seg2qc_note[segment_name][var] = seg_qc_note

          # compile segment qc notes to qc note
          t_qc_note = ""
          for var in segment_check_criteria.keys():
            for segment in seg2qc_note:
              t_qc_note += seg2qc_note[segment][var]
          if t_qc_note:
            qc_status = "QC_ALERT"
            qc_note += t_qc_note
          qc_check_metrics = [m for m in qc_check_metrics if m not in to_rm]

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
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}