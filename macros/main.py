import csv
import os
import fnmatch

### THIS IS A MACROS OVERRIDE FOR THE DOCUMENTATION ###


def define_env(env):
  @env.macro
  def input_table(filename=None, filter_column=None, filter_values=None, columns=None, sort_by=None, input_table=False, indent=0):
    if not filename:
      return "**Error:** `filename` is required."
    
    # Normalize filter_values list (if provided as a string, split by commas)
    if filter_values:
      if isinstance(filter_values, str):
          filter_values = [v.strip() for v in filter_values.split(',')]
    else:
      filter_values = None  # Disable filtering if not provided

    # Load CSV
    with open(filename, newline='') as csvfile:
      reader = csv.DictReader(csvfile, delimiter='\t')
      all_headers = reader.fieldnames

      # Check if the filter_column exists in the headers
      filter_column_exists = filter_column in all_headers if filter_column else False

      # Expand wildcard column patterns
      if columns:
        matched_columns = []
        for col in columns:
          matches = fnmatch.filter(all_headers, col)
          matched_columns.extend(matches)
        # Remove duplicates, preserve order
        seen = set()
        headers = [h for h in matched_columns if not (h in seen or seen.add(h))]
      else:
        headers = all_headers

      # Filter rows by filter_column (if specified and filter_column exists)
      rows = []
      for row in reader:
        if filter_column_exists:
          row_values = [v.strip() for v in row[filter_column].split(',')]
          if not filter_values or any(fv in row_values for fv in filter_values):
            rows.append(row)
        else:
          # If no filter column, include all rows
          rows.append(row)

    # Sort rows if requested
    if sort_by:
      if isinstance(sort_by, str):
        sort_by = [(sort_by, False)] # default to ascending sort order
        
      elif isinstance(sort_by, list):
        # Normalize: convert ["col1", "col2"] → [("col1", False), ("col2", False)]
        sort_by = [(col, False) if isinstance(col, str) else col for col in sort_by]

      for col, reverse in reversed(sort_by):
          rows.sort(key=lambda r: r.get(col, ''), reverse=reverse)


    indent_str = ' ' * indent
    # Build Markdown table
    md = '| ' + ' | '.join(f'**{h}**' for h in headers) + ' |\n'
    md += indent_str + '| ' + ' | '.join(['---'] * len(headers)) + ' |\n'
    # Iterate through the rows and bold the second column (index 1) if input table is True
    for row in rows:
      if input_table:
        md += indent_str + '| ' + ' | '.join(f'**{row.get(h, "")}**' if index == 1 else row.get(h, '') for index, h in enumerate(headers)) + ' |\n'
      else:
        md += indent_str + '| ' + ' | '.join(row.get(h, '') for h in headers) + ' |\n'
    
    return md
