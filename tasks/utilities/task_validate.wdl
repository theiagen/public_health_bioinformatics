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
    File datatable1_tsv
    File datatable2_tsv
    File validation_criteria_tsv
    String columns_to_compare

    Int disk_size = 10
  }
  command <<<

    python3 <<CODE
  import pandas as pd
  import numpy as np
  import os
  import sys
  import pdfkit as pdf

  def read_tsv(tsv_file):
    # Read TSV and change first column to 'samples'
    df = pd.read_csv(tsv_file, sep='\t')
    c1_name = df.columns.values[0]
    df.columns.values[0] = "samples"
    
    # replace blank cells with NaNs 
    df.replace(r'^\s+$', np.nan, regex=True)

    return [df, c1_name]

  # Read in TSVs and keep table name
  df1, df1_c1_name = read_tsv("~{datatable1_tsv}")
  df2, df2_c1_name = read_tsv("~{datatable2_tsv}")

  # read in the list of columns to compare
  comparison_columns = "~{columns_to_compare}".split(",")
  print(comparison_columns)

  # remove columns that we will not compare from the two data tables
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
  new columns_df2 = []

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
  
  # get count of nans per column
  df1_nan_counts = pd.DataFrame(df1.isna().sum(axis = 0), columns = ['sum'])
  df2_nan_counts = pd.DataFrame(df2.isna().sum(axis = 0), columns = ['sum'])
  

  # Perform comparison
  df_comp1 = df1.compare(df2, align_axis=1, keep_shape=True, keep_equal=False)
  # Count non-NA values in each columns
  val_cnts = df_comp1.count()
  df_val_cnts=val_cnts.to_frame()
  df_val_cnts.columns = ['Number of Diffs']
  print(df_val_cnts)
  # Replace NAs with "EXACT_MATCH"
  df_comp1.fillna(value='EXACT_MATCH', method=None, axis=None, inplace=True, limit=None, downcast=None)

  # Get the side-by-side comparison of the TSVs
  # df_diff_vert = df1.compare(df2, align_axis = 0, keep_shape=True, keep_equal=True).transpose()
  # df_comp_bool = df1.where()
  #missing_samples = []

  # Compare samples
  #print(df1.samples.values)
  #print('column\tis_same\ttsv1\ttsv2')
  #print(f'filename\t{args.tsv1==args.tsv2}\t{args.tsv1}\t{args.tsv2}')
  #for i, data in df_diff_vert.iterrows():
  #    if data[0]['self'] != data[0]['other']:
  #        print(f"{data.name}\t{data[0]['self'] == data[0]['other']}\t{data[0]['self']}\t{data[0]['other']}")


  count_dict={}
  for i in df_comp1.columns:
    count_dict[i]=df_comp1[i].value_counts()
  counts_df=pd.DataFrame.from_dict(count_dict, orient='columns', dtype=None, columns=None)
  print(counts_df)


  out_xlsx_name=f'{args.outdir}/{args.prefix}.xlsx'
  out_html_name=f'{args.outdir}/{args.prefix}.html'
  out_pdf_name=f'{args.outdir}/{args.prefix}.pdf'

  pd.set_option('display.max_colwidth', 20)
  df_comp1.to_excel(out_xlsx_name)
  df_val_cnts.to_html(out_html_name)

  options = {
  'page-width': '10000mm',
  'title': 'Validation Report',
  'margin-top': '0.25in',
  'margin-right': '0.25in',
  'margin-bottom': '0.25in',
  'margin-left': '0.25in'}
  out_pdf_var=out_html_name
  pdf1 = pdf.from_file(out_html_name, out_pdf_name, options=options)


  print(df_comp1)

  CODE
  >>>
  runtime {
    docker: "quay.io/theiagen/utility:1.2"
    memory: "4 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
  output {
    File pdf_report = "~{out_dir}/~{out_prefix}.pdf"
    File xl_report = "~{out_dir}/~{out_prefix}.xlsx"
  }
}