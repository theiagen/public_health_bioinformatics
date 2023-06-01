version 1.0

task export_two_tsvs {
  input {
    String terra_project
    String terra_workspace
    String datatable1
    String datatable2
    Int disk_size = 10
  }
  command <<<
    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable1} --tsv_filename ~{datatable1}

    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project} --workspace ~{terra_workspace} --entity_type ~{datatable2} --tsv_filename ~{datatable2}
  >>>
  runtime {
    docker: "quay.io/theiagen/terra-tools:2023-03-16"
    memory: "1 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File datatable1_tsv = "~{datatable1}"
    File datatable2_tsv = "~{datatable2}"
  }
}

task compare_two_tsvs {
  input {
    String datatable1
    File datatable1_tsv
    String datatable2
    File datatable2_tsv
    File validation_criteria_tsv
    String columns_to_compare
    String output_prefix

    Int disk_size = 10
  }
  command <<<

    pip install pretty_html_table

    python3 <<CODE
  import pandas as pd
  import numpy as np
  import os
  import sys
  import re
  import pdfkit as pdf
  from pretty_html_table import build_table

  def read_tsv(tsv_file):
    # Read TSV and change first column to 'samples'
    df = pd.read_csv(tsv_file, sep='\t')
    c1_name = df.columns.values[0]
    df.columns.values[0] = "samples"
    
    # replace blank cells with NaNs 
    df = df.replace(r'^\s+$', np.nan, regex=True)

    return [df, c1_name]

  # Read in TSVs and keep table name
  df1, df1_c1_name = read_tsv("~{datatable1_tsv}")
  df2, df2_c1_name = read_tsv("~{datatable2_tsv}")

  # read in the list of columns to compare
  comparison_columns = "~{columns_to_compare}".split(",")

  # add the "samples" column to the comparison_columns list
  comparison_columns.append("samples")

  # remove columns that we will not compare from the two data tables
  # initalize drop lists:
  drop_list1 = []
  drop_list2 = []

  for item in df1.columns:
    if item not in comparison_columns:
      drop_list1.append(item)
  df1.drop(drop_list1, axis='columns', inplace=True)

  for item in df2.columns:
    if item not in comparison_columns:
      drop_list2.append(item)
  df2.drop(drop_list2, axis='columns', inplace=True)

  # initialize a list of new columns found
  new_columns_df1 = []
  new_columns_df2 = []

  # identify any columns that do not appear in a data table
  for column in comparison_columns:
    if column not in df1.columns:
      new_columns_df1.append(column)
      # add the column to data table
      df1[column] = np.nan
    if column not in df2.columns:
      new_columns_df2.append(column)
      # add the column to data table
      df2[column] = np.nan
  
  # get count of populated cells per column
  df1_populated_rows = pd.DataFrame(df1.count(), columns = ['Number of samples populated in ~{datatable1}'])
  df2_populated_rows = pd.DataFrame(df2.count(), columns = ['Number of samples populated in ~{datatable2}'])
  
  # remove the sample name rows from the summary_output table (should be identical, no point checking here)
  df1_populated_rows.drop("samples", axis=0, inplace=True)
  df2_populated_rows.drop("samples", axis=0, inplace=True)

  # make a union of the two tables along the rows 
  summary_output = pd.concat([df1_populated_rows, df2_populated_rows], join="outer", axis=1)

  # count the number of differences using exact match
  # temporarily make NaNs Null since NaN != NaN for the pd.DataFrame.eq() function
  number_of_differences = pd.DataFrame((~df1.fillna("NULL").eq(df2.fillna("NULL"))).sum(), columns = ['Number of differences (exact match)'])
  # remove the sample name row 
  number_of_differences.drop("samples", axis=0, inplace=True)
  
  # add the number of differences to the summary output table
  summary_output = pd.concat([summary_output, number_of_differences], join="outer", axis=1)

  # get table of self-other differences
  # also: temporarily drop the sample name column for the comparison and then set it as the index for the output data frame
  df_comparison = df1.drop('samples', axis=1).compare(df2.drop('samples', axis=1), keep_shape=True).set_index(df1['samples'])
  # rename the self and other with table names
  df_comparison = df_comparison.rename(columns={'self': '~{datatable1}', 'other': '~{datatable2}'}, level=-1)
  # replace matching values (NAs) with blanks
  df_comparison.fillna('')

  # perform validation from validation tsv
  validation_criteria = pd.read_csv("~{validation_criteria_tsv}", sep='\t', index_col=0)
  validation_criteria = validation_criteria.transpose()
  # correct dtypes - convert to numeric first, and then to strings
  validation_criteria = validation_criteria.apply(pd.to_numeric, errors='ignore').convert_dtypes()

  # calculate percent difference with mean
  def percent_difference(col1, col2):
    # |x-y|/((x+y)/2)
    return np.absolute(col2.sub(col1)).div((col2.add(col1))/2)

  # perform validation checks
  def validate(series, df1, df2):
    if  pd.api.types.is_string_dtype(series) == True:
      if series[0] == "EXACT": # count number of exact match differences
        return ("EXACT", (~df1[series.name].fillna("NULL").eq(df2[series.name].fillna("NULL"))).sum())
      elif series[0] == "IGNORE": # do not check
        return ("IGNORE", 0)
      else:
        return("String value not recognized", np.nan)
    elif pd.api.types.is_float_dtype(series) == True:
      # calculate percent difference; compare percent difference to series[0] and return total count where it is greater than
      percent_difference(df1[series.name], df2[series.name]).gt(series[0]).sum()
      return(format(series[0], '.2%'), percent_difference(df1[series.name], df2[series.name]).gt(series[0]).sum())
    elif pd.api.types.is_datetime64_any_dtype(series) == True: # do not check
      return("DATE VALUE; IGNORED", 0)
    elif pd.api.types.is_integer_dtype(series) == True: # do something? we haven't decided yet.
      return("INTEGER; IGNORED FOR NOW", 0)
    else: # it's an object type
      return("OBJECT TYPE VALUE; IGNORED FOR NOW", 0)

  # perform check and add to the summary output table
  summary_output[["Validation Criteria", "Number of samples failing the validation criteria"]] = pd.DataFrame(validation_criteria.apply(lambda x: validate(x, df1, df2), result_type="expand")).transpose()

  print("successful 170")

  out_xlsx_name = "~{output_prefix}.xlsx"
  out_html_name = "~{output_prefix}.html"
  out_pdf_name = "~{output_prefix}.pdf"

  pd.set_option('display.max_colwidth', None)
  df_comparison.to_excel(out_xlsx_name)

  # make pretty html table
  html_table_light_grey = build_table(summary_output, 
    'grey_light', 
    index=True, 
    text_align='center',
    # conditions={
    #   'Number of differences (exact match)': {
    #     'min': 1,
    #     'max': 0,
    #     'min_color': 'black',
    #     'max_color': 'red'
    #   }
    # }
  )

  # save to html file
  with open(out_html_name, 'w') as outfile:
    outfile.write(html_table_light_grey)

  # convert to pdf
  options = {
    'page-size': 'Letter',
    'title': '~{datatable1} vs. ~{datatable2}',
    'margin-top': '0.25in',
    'margin-right': '0.25in',
    'margin-bottom': '0.25in',
    'margin-left': '0.25in'
  }

  output_pdf = pdf.from_file(out_html_name, out_pdf_name, options=options)

  CODE
  >>>
  runtime {
    docker: "quay.io/theiagen/utility:1.2"
    memory: "4 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File pdf_report = "~{output_prefix}.pdf"
    File excel_report = "~{output_prefix}.xlsx"
  }
}