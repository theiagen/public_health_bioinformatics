version 1.0

task qc_check {
  input {
    File? qc_check_table
    String? expected_taxon
    String? gambit_predicted_taxon
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
    String? midas_secondary_genus_abundance 
    Int? assembly_length
    Int? number_contigs 
    Int? n50_value 
    Float? ani_highest_percent 
    Float? ani_highest_percent_bases_aligned
    String? ani_top_species_match
    String? busco_results
    Int disk_size = 100
  }
  command <<<
    # date 
    date | tee DATE

    python3 <<CODE
    import csv
    import pandas as pd
    import numpy as np

    # create two empty variables that will be combined to create "qc_check" output string
    qc_status = ""
    qc_note = ""

    # read qc_check_table
    qc_check_df = pd.read_csv('~{qc_check_table}', sep='\t', index_col = 'taxon')
    qc_check_taxa = qc_check_df.index.values.tolist()

    # check that either expected_taxon or gambit_predicted_taxon is provided
    # use expected_taxon for qc_check, else use gambit_predicted_taxon
    expected_taxon = "~{expected_taxon}"
    gambit_predicted_taxon = "~{gambit_predicted_taxon}"
    if (expected_taxon):
      qc_taxon = expected_taxon.replace(" ", "_")
      print("User-provided expected_taxon was found, this taxon will be used for qc_check")
    elif (gambit_predicted_taxon):
      qc_taxon = gambit_predicted_taxon.replace(" ", "_")
      print("No user-provided expected_taxon was found, gambit_predicted_taxon will be used for qc_check")
    else:
      qc_status = str("QC_NA")
      qc_note = str(" No expected_taxon or gambit_predicted_taxon found, qc_check task could not proceed")

    # if qc_taxon variable was provided, check if qc_taxon is in qc_check_table
    if (qc_status != "QC_NA"):
      match = []
      for taxa in qc_check_taxa:
        if taxa in qc_taxon:
          match.append(taxa)
      if len(match) < 1: 
        qc_status = str("QC_NA")
        qc_note = str(" No matching taxon detected in qc_check_table")
      elif len(match) > 1:
        qc_status = str("QC_NA")
        qc_note = str(" Multiple matching taxa detected in qc_check_table")
      else: 
        taxon_df = qc_check_df.loc[[match[0]]]
        print(f"Exactly one matching taxon in qc_check table: {match}, proceeding with qc_check")

        # remove columns where all values are null
        taxon_df = taxon_df.replace(r'^\s*$', np.nan, regex=True)
        taxon_df = taxon_df.dropna(how='all', axis=1)
  
    ### perform QC checks for any metrics listed in qc_check_table
    if (qc_status != "QC_NA"):
      qc_check_metrics = taxon_df.columns.values.tolist()

      # check that est_coverage_raw is greater than metric
      est_coverage_raw = "~{est_coverage_raw}"
      if ("est_coverage_raw" in qc_check_metrics):
        if (est_coverage_raw):
          est_coverage_raw_metric = taxon_df['est_coverage_raw'][0]
          if (float(est_coverage_raw) >= est_coverage_raw_metric):
            print("est_coverage_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f" est_coverage_raw ({est_coverage_raw}) was less than minimum threshold of {est_coverage_raw_metric};"
          qc_check_metrics.remove("est_coverage_raw")
      else:
        print("est_coverage_raw not detected in qc_check_table")

      # check that est_coverage_clean is greater than metric
      est_coverage_clean = "~{est_coverage_clean}"
      if ("est_coverage_clean" in qc_check_metrics):
        if (est_coverage_clean):
          est_coverage_clean_metric = taxon_df['est_coverage_clean'][0]
          if (float(est_coverage_clean) >= est_coverage_clean_metric):
            print("est_coverage_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f" est_coverage_clean ({est_coverage_clean}) was less than minimum threshold of {est_coverage_clean_metric};"
          qc_check_metrics.remove("est_coverage_clean")
      else:
        print("est_coverage_clean not detected in qc_check_table")

      # check that r1_mean_q_raw is greater than metric
      r1_mean_q_raw = "~{r1_mean_q_raw}"
      if ("r1_mean_q_raw" in qc_check_metrics):
        if (r1_mean_q_raw):
          r1_mean_q_raw_metric = taxon_df['r1_mean_q_raw'][0]
          if (float(r1_mean_q_raw) >= r1_mean_q_raw_metric):
            print("r1_mean_q_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r1_mean_q_raw ({r1_mean_q_raw}) was less than minimum threshold of {r1_mean_q_raw_metric};"
          qc_check_metrics.remove("r1_mean_q_raw")
      else:
        print("r1_mean_q_raw not detected in qc_check_table")

      # check that r2_mean_q_raw is greater than metric
      r2_mean_q_raw = "~{r2_mean_q_raw}"
      if ("r2_mean_q_raw" in qc_check_metrics):
        if (r2_mean_q_raw):
          r2_mean_q_raw_metric = taxon_df['r2_mean_q_raw'][0]
          if (float(r2_mean_q_raw) >= r2_mean_q_raw_metric):
            print("r2_mean_q_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r2_mean_q_raw ({r2_mean_q_raw}) was less than minimum threshold of {r2_mean_q_raw_metric};"
          qc_check_metrics.remove("r2_mean_q_raw")
      else:
        print("r2_mean_q_raw not detected in qc_check_table")

      # check that combined_mean_q_raw is greater than metric
      combined_mean_q_raw = "~{combined_mean_q_raw}"
      if ("combined_mean_q_raw" in qc_check_metrics):
        if (combined_mean_q_raw):
          combined_mean_q_raw_metric = taxon_df['combined_mean_q_raw'][0]
          if (float(combined_mean_q_raw) >= combined_mean_q_raw_metric):
            print("combined_mean_q_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} combined_mean_q_raw ({combined_mean_q_raw}) was less than minimum threshold of {combined_mean_q_raw_metric};"
          qc_check_metrics.remove("combined_mean_q_raw")
      else:
        print("combined_mean_q_raw not detected in qc_check_table")        

      # check that r1_mean_readlength_raw is greater than metric
      r1_mean_readlength_raw = "~{r1_mean_readlength_raw}"
      if ("r1_mean_readlength_raw" in qc_check_metrics):
        if (r1_mean_readlength_raw):
          r1_mean_readlength_raw_metric = taxon_df['r1_mean_readlength_raw'][0]
          if (float(r1_mean_readlength_raw) >= r1_mean_readlength_raw_metric):
            print("r1_mean_readlength_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r1_mean_readlength_raw ({r1_mean_readlength_raw}) was less than minimum threshold of {r1_mean_readlength_raw_metric};"
          qc_check_metrics.remove("r1_mean_readlength_raw")
      else:
        print("r1_mean_readlength_raw not detected in qc_check_table")

      # check that r2_mean_readlength_raw is greater than metric
      r2_mean_readlength_raw = "~{r2_mean_readlength_raw}"
      if ("r2_mean_readlength_raw" in qc_check_metrics):
        if (r2_mean_readlength_raw):
          r2_mean_readlength_raw_metric = taxon_df['r2_mean_readlength_raw'][0]
          if (float(r2_mean_readlength_raw) >= r2_mean_readlength_raw_metric):
            print("r2_mean_readlength_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r2_mean_readlength_raw ({r2_mean_readlength_raw}) was less than minimum threshold of {r2_mean_readlength_raw_metric};"
          qc_check_metrics.remove("r2_mean_readlength_raw")
      else:
        print("r2_mean_readlength_raw not detected in qc_check_table")

      # check that combined_mean_readlength_raw is greater than metric
      combined_mean_readlength_raw = "~{combined_mean_readlength_raw}"
      if ("combined_mean_readlength_raw" in qc_check_metrics):
        if (combined_mean_readlength_raw):
          combined_mean_readlength_raw_metric = taxon_df['combined_mean_readlength_raw'][0]
          if (float(combined_mean_readlength_raw) >= combined_mean_readlength_raw_metric):
            print("combined_mean_readlength_raw passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} combined_mean_readlength_raw ({combined_mean_readlength_raw}) was less than minimum threshold of {combined_mean_readlength_raw_metric};"
          qc_check_metrics.remove("combined_mean_readlength_raw")
      else:
        print("combined_mean_readlength_raw not detected in qc_check_table")  

      # check that r1_mean_q_clean is greater than metric
      r1_mean_q_clean = "~{r1_mean_q_clean}"
      if ("r1_mean_q_clean" in qc_check_metrics):
        if (r1_mean_q_clean):
          r1_mean_q_clean_metric = taxon_df['r1_mean_q_clean'][0]
          if (float(r1_mean_q_clean) >= r1_mean_q_clean_metric):
            print("r1_mean_q_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r1_mean_q_clean ({r1_mean_q_clean}) was less than minimum threshold of {r1_mean_q_clean_metric};"
          qc_check_metrics.remove("r1_mean_q_clean")
      else:
        print("r1_mean_q_clean not detected in qc_check_table")

      # check that r2_mean_q_clean is greater than metric
      r2_mean_q_clean = "~{r2_mean_q_clean}"
      if ("r2_mean_q_clean" in qc_check_metrics):
        if (r2_mean_q_clean):
          r2_mean_q_clean_metric = taxon_df['r2_mean_q_clean'][0]
          if (float(r2_mean_q_clean) >= r2_mean_q_clean_metric):
            print("r2_mean_q_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r2_mean_q_clean ({r2_mean_q_clean}) was less than minimum threshold of {r2_mean_q_clean_metric};"
          qc_check_metrics.remove("r2_mean_q_clean")
      else:
        print("r2_mean_q_clean not detected in qc_check_table")

      # check that combined_mean_q_clean is greater than metric
      combined_mean_q_clean = "~{combined_mean_q_clean}"
      if ("combined_mean_q_clean" in qc_check_metrics):
        if (combined_mean_q_clean):
          combined_mean_q_clean_metric = taxon_df['combined_mean_q_clean'][0]
          if (float(combined_mean_q_clean) >= combined_mean_q_clean_metric):
            print("combined_mean_q_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} combined_mean_q_clean ({combined_mean_q_clean}) was less than minimum threshold of {combined_mean_q_clean_metric};"
          qc_check_metrics.remove("combined_mean_q_clean")
      else:
        print("combined_mean_q_clean not detected in qc_check_table")        

      # check that r1_mean_readlength_clean is greater than metric
      r1_mean_readlength_clean = "~{r1_mean_readlength_clean}"
      if ("r1_mean_readlength_clean" in qc_check_metrics):
        if (r1_mean_readlength_clean):
          r1_mean_readlength_clean_metric = taxon_df['r1_mean_readlength_clean'][0]
          if (float(r1_mean_readlength_clean) >= r1_mean_readlength_clean_metric):
            print("r1_mean_readlength_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r1_mean_readlength_clean ({r1_mean_readlength_clean}) was less than minimum threshold of {r1_mean_readlength_clean_metric};"
          qc_check_metrics.remove("r1_mean_readlength_clean")
      else:
        print("r1_mean_readlength_clean not detected in qc_check_table")

      # check that r2_mean_readlength_clean is greater than metric
      r2_mean_readlength_clean = "~{r2_mean_readlength_clean}"
      if ("r2_mean_readlength_clean" in qc_check_metrics):
        if (r2_mean_readlength_clean):
          r2_mean_readlength_clean_metric = taxon_df['r2_mean_readlength_clean'][0]
          if (float(r2_mean_readlength_clean) >= r2_mean_readlength_clean_metric):
            print("r2_mean_readlength_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} r2_mean_readlength_clean ({r2_mean_readlength_clean}) was less than minimum threshold of {r2_mean_readlength_clean_metric};"
          qc_check_metrics.remove("r2_mean_readlength_clean")
      else:
        print("r2_mean_readlength_clean not detected in qc_check_table")

      # check that combined_mean_readlength_clean is greater than metric
      combined_mean_readlength_clean = "~{combined_mean_readlength_clean}"
      if ("combined_mean_readlength_clean" in qc_check_metrics):
        if (combined_mean_readlength_clean):
          combined_mean_readlength_clean_metric = taxon_df['combined_mean_readlength_clean'][0]
          if (float(combined_mean_readlength_clean) >= combined_mean_readlength_clean_metric):
            print("combined_mean_readlength_clean passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} combined_mean_readlength_clean ({combined_mean_readlength_clean}) was less than minimum threshold of {combined_mean_readlength_clean_metric};"
          qc_check_metrics.remove("combined_mean_readlength_clean")
      else:
        print("combined_mean_readlength_clean not detected in qc_check_table")  

      # check that midas_secondary_genus_abundance is lower than metric
      midas_secondary_genus_abundance = "~{midas_secondary_genus_abundance}"
      if ("midas_secondary_genus_abundance" in qc_check_metrics):
        if (midas_secondary_genus_abundance):
          midas_secondary_genus_abundance_metric = taxon_df['midas_secondary_genus_abundance'][0]
          if (float(midas_secondary_genus_abundance) < midas_secondary_genus_abundance_metric):
            print("midas_secondary_genus_abundance passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} midas_secondary_genus_abundance ({midas_secondary_genus_abundance}) was greater than maximum threshold of {midas_secondary_genus_abundance_metric};"
          qc_check_metrics.remove("midas_secondary_genus_abundance")
      else:
        print("midas_secondary_genus_abundance not detected in qc_check_table")

      # check that assembly length is within acceptable range
      assembly_length = "~{assembly_length}"
      if ("assembly_length_min" in qc_check_metrics) and ("assembly_length_max" in qc_check_metrics):
        if (assembly_length):
          qc_check_metrics.remove("assembly_length_min")
          qc_check_metrics.remove("assembly_length_max")
          assembly_length_min_metric = taxon_df['assembly_length_min'][0]
          assembly_length_max_metric = taxon_df['assembly_length_max'][0]
          if (int(assembly_length) >= assembly_length_min_metric) and (int(assembly_length) <= assembly_length_max_metric):
            print("assembly_length passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} assembly_length ({assembly_length}) was outside of acceptable range ({assembly_length_min_metric} to {assembly_length_max_metric});"
      else:
        print("Either assembly_length_min or assembly_length_max was not detected in qc_check_table")

      # check that number_contigs is less than metric
      number_contigs = "~{number_contigs}"
      if ("number_contigs" in qc_check_metrics):
        if (number_contigs):
          number_contigs_metric = taxon_df['number_contigs'][0]
          if (float(number_contigs) <= number_contigs_metric):
            print("number_contigs passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} number_contigs ({mumber_contigs}) was greater than maximum threshold of {number_contigs_metric};"
          qc_check_metrics.remove("number_contigs")
      else:
        print("number_contigs not detected in qc_check_table")

      # check that n50_value is greater than metric
      n50_value = "~{n50_value}"
      if ("n50_value" in qc_check_metrics):
        if (n50_value):
          n50_value_metric = taxon_df['n50_value'][0]
          if (float(n50_value) >= n50_value_metric):
            print("n50_value passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} n50_value ({n50_value}) was less than minimum threshold of {n50_value_metric};"
          qc_check_metrics.remove("n50_value")
      else:
        print("n50_value not detected in qc_check_table")
        
      # check that ani_highest_percent is greater than metric
      ani_highest_percent = "~{ani_highest_percent}"
      if ("ani_highest_percent" in qc_check_metrics):
        if (ani_highest_percent):
          ani_highest_percent_metric = taxon_df['ani_highest_percent'][0]
          if (float(ani_highest_percent) >= ani_highest_percent_metric):
            print("ani_highest_percent passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} ani_highest_percent ({ani_highest_percent}) was less than minimum threshold of {ani_highest_percent_metric};"
          qc_check_metrics.remove("ani_highest_percent")
      else:
        print("ani_highest_percent not detected in qc_check_table")

      # check that ani_highest_percent_bases_aligned is greater than metric
      ani_highest_percent_bases_aligned = "~{ani_highest_percent_bases_aligned}"
      if ("ani_highest_percent_bases_aligned" in qc_check_metrics):
        if (ani_highest_percent_bases_aligned):
          ani_highest_percent_bases_aligned_metric = taxon_df['ani_highest_percent_bases_aligned'][0]
          if (float(ani_highest_percent_bases_aligned) >= ani_highest_percent_bases_aligned_metric):
            print("ani_highest_percent_bases_aligned passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} ani_highest_percent_bases_aligned ({ani_highest_percent_bases_aligned}) was less than minimum threshold of {ani_highest_percent_bases_aligned_metric};"
          qc_check_metrics.remove("ani_highest_percent_bases_aligned")
      else:
        print("ani_highest_percent_bases_aligned not detected in qc_check_table")

      # check that busco_results is greater than metric
      busco_results = "~{busco_results}"
      if ("busco_completeness" in qc_check_metrics):
        if (busco_results):
          busco_completeness = busco_results.split("C:")[1].split("%")[0]
          qc_check_metrics.remove("busco_completeness")
          busco_completeness_metric = taxon_df['busco_completeness'][0]
          if (float(busco_completeness) >= busco_completeness_metric):
            print("busco_completeness passed qc_check")
          else:
            qc_status = "QC_ALERT"
            qc_note = f"{qc_note} busco_completeness ({busco_completeness}) was less than minimum threshold of {busco_completeness_metric};" 
      else:
        print("busco_completeness not detected in qc_check_table")

      # if any columns in the qc_check_table went unused, issue a QC ALERT
      if (len(qc_check_metrics) > 0):
        qc_status = "QC_ALERT"
        qc_note = f"{qc_note} one or more columns in qc_check_table was missing a required input: {qc_check_metrics}"

    # if there was no QC_NA or QC_ALERT, then the sample is considered a QC_PASS
    if (qc_status != "QC_NA") and (qc_status != "QC_ALERT"):
      qc_check = "QC_PASS"
    # otherwise, report status and reason
    else:
      qc_check = str(qc_status) + ":" + str(qc_note)
      qc_check = qc_check.rstrip(';')
    with open("QC_CHECK", 'wt') as f:
        f.write(str(qc_check))

    CODE

  >>>
  output {
    String qc_check = read_string("QC_CHECK")
    File? qc_standard = qc_check_table
    String date = read_string("DATE")
  }
  runtime {
    docker: "broadinstitute/terra-tools:tqdm"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}