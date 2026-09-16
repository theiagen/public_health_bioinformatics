import csv
import os
import re
import fnmatch
import ast
from pathlib import Path
from functools import lru_cache

# cache files to prevent rereading
@lru_cache(maxsize=256)
def read_file(path):
    return Path(path).read_text(encoding="utf-8")

##### RENDER TSV TABLE #####

# parse filtering parameters for rendering tsv
def normalize_filters(filters):    
    """
    Converts filter input into a consistent format:
    { column_name: [list of values] }

    This ensures downstream filtering logic always works with lists,
    regardless of whether the input was a string, list, or single value.
    """
    if not filters:
        return {}
    out = {}

    # Case 1: filter value is a string; Example: "A,B,C"
    # Case 2: filter is list or tuple; Example: ["A", "B", "C"]
    # Case 3: single scalar value (int, bool, etc.); Example: 5 or True
    for col, val in filters.items():
        if isinstance(val, str):
            # Split by comma and remove extra whitespace; discard empty values
            out[col] = [v.strip() for v in val.split(",") if v.strip()]
        elif isinstance(val, (list, tuple)):
            # Convert all values to strings and strip whitespace
            out[col] = [str(v).strip() for v in val]
        else:
            out[col] = [str(val)]
    return out

# extract rows that match filters
def row_matches(row, filters):    
    """
    Checks whether a single data row satisfies all filter conditions.

    A row is included only if:
    - Every filter column exists in the row
    - At least one allowed value matches the row's cell values
    """
    # allowed = list of acceptable values for that column
    for col, allowed in filters.items():
        if col not in row:
            return False

        # Split the cell value into a list in case it contains multiple
        # comma-separated entries (e.g. "A,B,C")
        cell = [v.strip() for v in row[col].split(",") if v.strip()]

        if not any(v in cell for v in allowed):
            return False
    return True

# ensure sorting configuration works
def normalize_sort(sort_by):
    """
    Normalizes sorting configuration into a consistent format:
    [(column_name, reverse_flag), ...]

    This ensures downstream sorting logic always works with a list of
    (column, ascending/descending) tuples.
    """
    if not sort_by:
        return []

    if isinstance(sort_by, str): # single string; ascending sort by default
        return [(sort_by, False)]

    result = []

    for item in sort_by: # sort instructions
        if isinstance(item, str):  # column name; ascending
            result.append((item, False))
        elif isinstance(item, (tuple, list)) and len(item) == 2: #tuple with reverse flag
            result.append((item[0], bool(item[1])))
        else:
            raise ValueError(f"bad sort_by value: {item}")

    return result

# add tsv table to markdown
def render_tsv_table(
        filename=None, # name of table
        filters=None, # keep matching rows
        exclude_filters=None, # remove matching rows
        columns=None, # columns to include; wildcards allowed
        sort_by=None, # sorting order
        input_table=False, # apply specific formatting
        indent=0, # add indentations
        replacements=None, # replace certain language
    ):
    """
    Renders a TSV file into a Markdown table with:
    - column filtering
    - row filtering (include/exclude rules)
    - sorting
    - optional formatting tweaks
    """
    if not filename:
        return "**Error:** filename required"

    with open(filename, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    headers = list(rows[0].keys()) if rows else []

    # keep only desired columns
    if columns:
        picked = []
        for c in columns:
            picked.extend(fnmatch.filter(headers, c))
        headers = list(dict.fromkeys(picked))

    # apply filters and keep them in matching
    filters = normalize_filters(filters)
    rows = [r for r in rows if row_matches(r, filters)]

    # apply exclude filters and remove them if matching
    exclude_filters = normalize_filters(exclude_filters)
    if exclude_filters:
        # returns true if row should be excluded
        def should_remove(row): 
            for col, vals in exclude_filters.items():
                if col in row:
                    cell = [v.strip() for v in row[col].split(",")]
                    if all(v in cell for v in vals):
                        return True
            return False

        rows = [r for r in rows if not should_remove(r)]

    # sort the columns
    for col, rev in reversed(normalize_sort(sort_by)):
        rows.sort(key=lambda r: r.get(col, ""), reverse=rev)

    # create indentation string 
    indent_str = " " * indent

    md = "| " + " | ".join(headers) + " |\n"
    md += indent_str + "| " + " | ".join(["---"] * len(headers)) + " |\n"

    for row in rows:
        vals = []
        for index, header in enumerate(headers):
            value = row.get(header, "")
            if input_table and index == 1:
                value = f"**{value}**"

            vals.append(value)

        md += indent_str + "| " + " | ".join(vals) + " |\n"

    if replacements:
        for old, new in replacements.items():
            md = md.replace(old, new)

    return md

############################

######## INCLUDE MD ########

# parse include_md(...) arguments using AST (abstract syntax tree)
def parse_macro_args(argument_string):
    """
    Takes the include_md arguments (`common_text/file.md, indent=4, condition="theiaprok"`)
    and turns it into a dictionary:
    { "path": "common_text/file.md", "indent": 4, "condition": "theiaprok" }
    """
    try:
        # wrap raw macro args into a synthetic Python expression
        # so AST parsing can safely extract structure
        ast_wrapped_input = f"f({argument_string})"

        # ast.parse converts this into a tree structure
        parsed_ast = ast.parse(ast_wrapped_input, mode="eval")
        call = parsed_ast.body
        if not isinstance(call, ast.Call):
            raise ValueError("the argument string has bad macro syntax")

        parsed_args = {}

        if call.args:
            parsed_args["path"] = ast.literal_eval(call.args[0])

        for keyword in call.keywords:
            parsed_args[keyword.arg] = ast.literal_eval(keyword.value)

        return parsed_args

    except Exception as e:
        raise ValueError(f"couldn't parse macro args: {e}")

"""
Matches conditional "if" blocks inside Markdown comments:
  <!-- if: theiaprok|theiaviral -->

Captures the allowed condition names (split by "|") so that
the rendering engine can decide whether to include or skip
the following block of content.
"""
IF_RE = re.compile(r'<!--\s*if:\s*([\w|]+)\s*-->')

""" 
Matches the closing tag for a conditional block in Markdown:
  <!-- endif -->

Used to signal the end of a conditionally included section,
allowing the parser to restore normal rendering behavior.
"""
ENDIF_RE = re.compile(r'<!--\s*endif\s*-->')

class ConditionStack:
    """
    Tracks whether a line of Markdown should be included or excluded
    based on conditional blocks like:

        <!-- if: conditionA|conditionB -->
        ...
        <!-- endif -->

    The stack allows nested conditional blocks to be evaluated correctly.
    """
    def __init__(self, active):
        # The currently active workflow/condition (e.g. "theiaprok")
        # This is compared against allowed conditions in IF blocks.
        self.active = active
        
        # Stack representing nested conditional states.
        # Each entry is a boolean:
        #   True  -> content in this block is included
        #   False -> content in this block is skipped
        #
        # Start with True because outside any IF block,
        # content is always included by default.
        self.stack = [True]

    def update(self, line):   
        """
        Processes a single line of Markdown and updates conditional state.

        Returns:
            - True  → line should be included
            - False → line should be excluded
            - None  → line is a control directive (if/endif)
        """
        match = IF_RE.match(line)
        if match:
            allowed_conditions = [condition.strip() for condition in match.group(1).split("|")]
            is_included = self.active in allowed_conditions 
            
            self.stack.append(is_included)
            return None

        # Detect end of conditional block: <!-- endif -->
        if ENDIF_RE.match(line):
            # Only pop if we are inside a nested block
            # (prevents popping the base True state)
            if len(self.stack) > 1:
                self.stack.pop()
            return None

        # For normal content lines:
        # Only include if ALL active conditions are True
        return all(self.stack)

# adjust markdown headings and add indentations
def adjust_heading(line, offset, indent_str):    
    """
    Adjusts Markdown heading levels and optionally indents the line.

    Example:
        "# Title" with offset=1 → "## Title"
    """
    # Try to detect a Markdown heading:
    # ^            → start of line
    # (#{1,6})     → 1 to 6 '#' characters (heading level)
    # \s+          → at least one space
    # (.*)         → the rest of the text (heading title)
    match = re.match(r'^(#{1,6})\s+(.*)', line)
    # if not a match, just add the identation
    if not match:
        return indent_str + line

    # If it IS a heading, split it into parts:
    # hashes = "###"
    # text   = "My Heading"
    hashes, text = match.groups()
    new_level = len(hashes) + (offset or 0)
    # keep heading level between 1 and 6 (markdown only supports h1–h6)
    new_level = max(1, min(new_level, 6))

    # add indentation and new heading
    return indent_str + ("#" * new_level) + " " + text

# resolve relative links
def resolve_links(line, current_file, root_file):    
    """
    Rewrites relative Markdown links so they remain valid
    after files are included into a larger documentation structure.
    """
    def replace(match):
        is_img, text, url = match.groups()

        # leave external stuff alone
        if url.startswith(("http://", "https://", "#", "mailto:")):
            return match.group(0)

        # convert to absolute filepath
        target = (current_file.parent / url).resolve()

        # make relative to docs root
        try:
            rel = os.path.relpath(target, root_file.parent)
        except Exception:
            # if something breaks just keep original
            return match.group(0)

        rel = rel.replace(os.sep, "/")

        # keep the ! if it was an image
        prefix = "!" if is_img else ""
        return f"{prefix}[{text}]({rel})"

    # add replacement to all markdown links in the line 
    # looks for [text](link) syntax with optional ! at front
    return re.sub(r'(!)?\[(.*?)\]\((.*?)\)', replace, line)


# match include md jinja calls
INCLUDE_RE = re.compile(r'{{\s*include_md\((.*?)\)\s*}}')

def include_md(
        path, # path to md file
        level_offset=None, # offset headings by x amount
        indent=0, # add x indentation
        condition=None, # content to include conditionals
        replacements=None, # text replacement
        _visited=None,
        _root=None,
    ):   
    """
    Recursively includes and processes Markdown files with support for:
    - nested includes
    - conditional blocks
    - heading adjustment
    - link resolution
    - indentation control
    - circular include protection
    """
    # resolve path relative to root
    docs = Path("docs").resolve()
    file_path = (docs / Path(path)).resolve()
    
    # create indentation string 
    indent_str = " " * indent

    # prevent infinite recursion
    if _visited is None:
        _visited = set()

    if file_path in _visited:
        raise RuntimeError(f"circular include: {file_path}")

    _visited.add(file_path)

    # make sure included file exists
    if not file_path.exists():
        return f"**Error:** missing file `{file_path}`"

    text = read_file(str(file_path))
    lines = text.splitlines()

    # remove YAML frontmatter (--- ... ---)
    if lines and lines[0].strip() == "---":
        for i in range(1, len(lines)):
            if lines[i].strip() == "---":
                lines = lines[i+1:]
                break

    # process conditions
    condition = ConditionStack(condition)
    out = []
    in_code_block = False

    for line in lines:
        # check conditional state
        state = condition.update(line)
        if state is None or not state:
            continue

        # check if this is a code block
        if line.lstrip().startswith("```"):
            in_code_block = not in_code_block
                
        # don't do any fancy with code blocks
        if not in_code_block:
            # see if we have a nested includes
            nested_match = INCLUDE_RE.match(line)
            if nested_match:
                # parse arguments inside macro call and inherit parent conditions
                args = parse_macro_args(nested_match.group(1))
                args.setdefault("condition", condition)
                # increase indentation relative to parent include
                args["indent"] = indent + args.get("indent", 0)

                # process recursively
                out.append(include_md(
                    _visited=_visited,
                    _root=_root,
                    **args
                ))
                continue

            if line.startswith("#"):
                out.append(adjust_heading(line, level_offset or 0, indent_str))
            else:
                # resolve links & only indent if line isn't empty
                resolved = resolve_links(line, file_path, _root or file_path)
                out.append(indent_str + resolved if resolved.strip() else "")
            
            continue
        
        # adjust within code block
        out.append(indent_str + line)
        
    # return output as single string
    result = "\n".join(out)

    # add post-processing text replacements
    if replacements:
        for old, new in replacements.items():
            result = result.replace(old, new)
    
    return result

# zensical entrypoint
def define_env(env):    
    """
    Registers custom macros into the Zensical environment.

    This function is called automatically by the macros plugin
    when the documentation build environment is initialized.
    """
    env.macro(include_md)
    env.macro(render_tsv_table)