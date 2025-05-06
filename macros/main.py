import csv
import os
import re
import fnmatch
import ast
from pathlib import Path
from urllib.parse import urljoin
from shlex import split as shlex_split

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
  
  @env.macro
  def include_md(path, level_offset=None, indent=0, condition=None, replacements=None):
    """Include a Markdown file, adjust heading levels and indent."""

    # Define the directory where the markdown files are located
    docs_dir = Path(env.project_dir) / 'docs'
    include_path = (docs_dir / Path(path)).resolve()

    # print(f"Resolved include_path: {include_path} — Exists? {include_path.is_file()}")
    # Check if the file exists, return an error if not
    if not include_path.exists():
      return f"**Error**: File not found: `{include_path}`"

    # Read the content of the included markdown file
    content = include_path.read_text(encoding='utf-8')
    lines = content.split('\n')
    adjusted_lines = []
    include_block = True

    # Iterate through the lines, adjusting heading levels and handling indentation
    for line in lines:
      # detect start of conditional block
      if_match = re.match(r'<!--\s*if:\s*([\w|]+)\s*-->', line)
      endif_match = re.match(r'<!--\s*endif\s*-->', line)

      if if_match:
        allowed = [x.strip() for x in if_match.group(1).split('|')]
        include_block = condition in allowed
        continue
      elif endif_match:
        include_block = True
        continue
      
      if include_block: 
        # Nested include detection
        nested_match = re.match(r'{{\s*include_md\((.*?)\)\s*}}', line)
        if nested_match:
          try:
            args_str = nested_match.group(1)
            args = parse_macro_args(args_str)
            args.setdefault("condition", condition)
            args["indent"] = indent + int(args.get('indent', 0))
            
            nested_content = include_md(**args)
            adjusted_lines.append(nested_content)      
          except Exception as e:
            adjusted_lines.append(f"**Error**:Nested include error: {e}")
            print(f"Error including nested file: {e}")
          continue

        if line.startswith('#'):
          adjusted_lines.append(adjust_heading(line, level_offset, indent))
        else:
          adjusted_lines.append(' ' * indent + resolve_links(line, include_path))

    result = '\n'.join(adjusted_lines)
  
    if replacements is not None:
      for old, new in replacements.items():
        result = result.replace(old, new)

    return result

  def adjust_heading(line, level_offset, indent):
    match = re.match(r'^(#{1,6})\s*(.*)', line)
    if match:
      hashes, text = match.groups()
      new_level = max(1, min(len(hashes) + (level_offset or 0), 6))
      adjusted_heading = ' ' * indent + new_level + ' ' + text
      return adjusted_heading
    else:
      return ' ' * indent + line
    
  def resolve_links(line, base_path):
    link_pattern = re.compile(r'(!)?\[(.*?)\]\((.*?)\)')  
    
    def replace_link(match):
      is_image = match.group(1) == '!'
      alt_text = match.group(2)
      url = match.group(3)

      # Only rewrite relative URLs (not http:// or https://)
      if url.startswith('http://') or url.startswith('https://'):
        return match.group(0)

      # Resolve the relative URL to the base path (original markdown file's location)
      absolute_path = (Path(base_path).parent / url).resolve()
      docs_dir = (Path(env.project_dir) / 'docs').resolve()
      
      try:
        resolved_path = os.path.relpath(absolute_path, start=Path(base_path).parent).replace(os.sep, '/')
      except ValueError:
        return match.group(0)
      
      # Return the correctly formatted link or image tag
      if is_image:
        return f'![{alt_text}]({resolved_path})'
      else:
        return f'[{alt_text}]({resolved_path})'

    return link_pattern.sub(replace_link, line)
    
  def parse_macro_args(arg_str):
    args = {}
    tokens = [token.strip() for token in arg_str.split(",")]
    for i, token in enumerate(tokens):      
      if token.startswith("replacements="):
        try:
          replacements_str = "=".join(token.split("=")[1:]).strip()
          replacements = ast.literal_eval(replacements_str)
          if not isinstance(replacements, dict):
            raise ValueError("Replacements must be a dictionary.")
          args["replacements"] = replacements
        except Exception as e:
          raise ValueError(f"Invalid replacements format: {e}")
    
      elif "=" in token:
        key, value = token.split("=", 1)
        val = value.strip().rstrip(',')
        
        if val.startswith('"') and val.endswith('"') or val.startswith("'") and val.endswith("'"):
          val = val[1:-1]
        
        elif val.isdigit():
          val = int(val)
        args[key.strip()] = val
      
      elif i == 0: 
        args['path'] = token.strip('"').strip("'").rstrip(',')
      
      else:
        raise ValueError(f"Invalid macro argument: {token}")
    return args