import csv
import os
import re
import fnmatch
from pathlib import Path
from urllib.parse import urljoin

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
  def include_md(path, level_offset=None, indent=0, condition=None):
    """Include a Markdown file, adjust heading levels and indent."""
    
    # Define the directory where the markdown files are located
    docs_dir = Path(env.project_dir) / 'docs'
    include_path = (docs_dir / path).resolve()

    # Check if the file exists, return an error if not
    if not include_path.exists():
      return f"**Error**: File not found: `{path}`"

    # Read the content of the included markdown file
    content = include_path.read_text(encoding='utf-8')

    # Split the content into lines for processing
    lines = content.split('\n')
    adjusted_lines = []

    # Iterate through the lines, adjusting heading levels and handling indentation
    for line in lines:
      # detect start of conditional block
      start_match = re.match(r'<!--\s*if:\s*(\w+)\s*-->', line)
      end_match = re.match(r'<!--\s*endif\s*-->', line)

      if start_match:
          if start_match.group(1) == condition:
              include = True
          else:
              include = False
          continue
      elif end_match:
          include = True
          continue     
        
      if include: 
        # Check if the line starts with a heading (i.e., begins with '#')
        if line.startswith('#'):
          # Adjust the heading by the specified level_offset and respect the indentation
          adjusted_lines.append(adjust_heading(line, level_offset, indent))
        else:
          # For non-heading lines, simply add indentation
          adjusted_lines.append(' ' * indent + resolve_links(line))

    # Join the adjusted lines into a single string and return the result
    return '\n'.join(adjusted_lines)

def adjust_heading(line, level_offset, indent):
  """Adjust heading level and maintain indentation."""
  # Match headings of the form '#', '##', '###', etc.
  match = re.match(r'^(#{1,6})\s*(.*)', line)
  
  if match:
    # Extract the current heading level (number of '#') and the heading text
    hashes = match.group(1)
    heading_text = match.group(2)
    
    # Calculate the new heading level, ensuring it's within valid bounds (1-6)
    new_level = max(1, min(len(hashes) + level_offset, 6))

    # Apply the indentation and adjust the heading level
    adjusted_heading = ' ' * indent + f"{'#' * new_level} {heading_text}"

    # Return the adjusted heading
    return adjusted_heading
  else:
    # Return the original line if it's not a heading
    return line
  
def resolve_links(line):
  """Resolve relative links in the markdown content."""
  # Regex to match Markdown links: [alt text](url)
  link_pattern = re.compile(r'!\[([^\]]*)\]\(([^)]+)\)|\[([^\]]+)\]\(([^)]+)\)')
  
  # Function to resolve the link based on the relative path
  def replace_link(match):
    # Extract the link components
    alt_text = match.group(1) or match.group(3)  # alt text for images or links
    url = match.group(2) or match.group(4)  # the URL in the markdown
    
    # Here, you can add logic to resolve relative paths (e.g., to your `docs/` folder)
    resolved_url = resolve_relative_path(url)

    # Return the updated markdown link
    return f"[{alt_text}]({resolved_url})" if alt_text else f"![{alt_text}]({resolved_url})"
  
  # Replace all links in the line
  return link_pattern.sub(replace_link, line)

def resolve_relative_path(url):
  """Resolve relative links to full paths."""
  # If the link is relative, we can adjust it based on the current document's location
  if url.startswith('http://') or url.startswith('https://'):
    return url  # No change for absolute URLs
  else:
    # For relative paths, resolve them relative to the current document's directory
    base_path = Path(url)
    return base_path.resolve()  # Use absolute path resolution for relative links
    