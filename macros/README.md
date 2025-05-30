# Documentation Macros

The macros in this module are used to:

- Render filtered and styled tables from `.tsv` (tab-separated) files
- Include and format Markdown content with heading adjustments, conditional blocks, and nested includes

These tools support modular and consistent documentation.

## Usage

### 1. `render_tsv_table(...)`

This macro renders a TSV file as a Markdown table. It supports column filtering, wildcard selection, sorting, and optional formatting for inputs.

#### `render_tsv_table()` Parameters

| Name            | Type                  | Required | Description |
|-----------------|-----------------------|----------|-------------|
| `filename`      | `str`                 | Yes      | Path to the TSV file (relative to the project root). |
| `filters`       | `Dict[str: str]`      | No       | Dictionary for fields to include (filters for, not out) {"Column name": "Target Value"} |
| `columns`       | `list[str]`           | No       | Subset of columns (supports wildcards). If not specified, all columns are shown. |
| `sort_by`       | `str`, `list[str]`, or `list[tuple]` | No | Sort by one or more columns. Use `(column, True)` for descending sort. |
| `input_table`   | `bool`                | No       | If `True`, bolds the second column (used for emphasizing input names). |
| `indent`        | `int`                 | No       | Number of spaces to indent the output Markdown. This will likely be a multiple of 4. |

#### `render_tsv_table()` Examples

Note that these examples include the `/// html | div[class="searchable-table"]` and `///` lines. These are required to enable the table-specific search bars.

##### Example 1: Input Tables

The following code is how this function should be called for creating the Inputs table.

The `columns` and `sort_by` parameters should be used for all input tables, and the `input_table` parameter should be set to `True` for all input tables.

```md
/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "TheiaProk_Illumina_SE"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"], indent=8) }}

///
```

##### Example 2: Output Tables

The following code is how this function should be called for creating the Outputs table.

The `columns` and `sort_by` parameters should be used for all output tables, and the `input_table` parameter should be set to `False` for all output tables.

```md
/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "TheiaProk_Illumina_PE"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"], indent=4) }}

///
```

---

### 2. `include_md(...)`

This macro includes and optionally transforms a Markdown file. It adjusts heading levels, resolves relative links, and supports conditional and nested includes.

#### `include_md()` Parameters

| Name           | Type                | Required | Description |
|----------------|---------------------|----------|-------------|
| `path`         | `str`               | Yes      | Path to the Markdown file to include (relative to `docs/`). |
| `level_offset` | `int`               | No       | Shifts heading levels (e.g., `#` → `##`) to match nesting context. UNTESTED, USE AT OWN RISK |
| `indent`       | `int`               | No       | Number of spaces to indent each line. |
| `condition`    | `str`               | No       | Include blocks only if they match this condition (e.g., `ONT`). |
| `replacements` | `dict[str, str]`    | No       | Dictionary of literal string replacements in the included content. |

#### `include_md()` Examples

##### Example 1: Basic Include

This example shows how to include a Markdown file  without any additional parameters.

```jinja
{{ include_md("common_text/bwa_task.md") }}
```

##### Example 2: Adding Indentation, a Condition, and Replacements Dictionary

This example shows how to add a replacements dictionary to replace the `??? task` string with `??? toggle` in the included file. This is because in the included documentation, the placement of this task required a toggle admonition instead of a task admonition, which is what is used in the individual task Markdown file.

```jinja
{{ include_md("common_text/kraken2_task.md", indent=4, condition="theiacov", replacements={"??? task": "??? toggle"}) }}
```

Note: Make sure any **assets** in the included Markdown file are prefixed with `../../` to ensure they are correctly resolved relative to the final output location. This is important for images, files, and othrr assets to ensure they are correctly displayed in the final documentation.

#### Conditional Block Syntax

You can include blocks conditionally using comment markers:

```markdown
<!-- if:ONT|SE -->
This content will appear only if condition is "ONT" or "SE".
<!-- endif -->
```

---

## For Analysts

1. Make sure to install the plugin with `pip install mkdocs-macros-plugin` before running `mkdocs serve`
2. Use `render_tsv_table()` to render TSV tables dynamically, avoiding copy-paste duplication. Add the TSV table to `docs/assets/tables/`. If you are adding inputs and outputs, make sure to add them to the `all_inputs.tsv` and `all_outputs.tsv` files, respectively. If you have an input and output that already exists for a different workflow, add your new workflow to the comma-separated list in the `Workflow` column of the TSV file. If you are adding a new input or output that does not exist in the other workflow, add it to the bottom of the appropriate TSV file. The `render_tsv_table()` function will automatically filter the table based on the workflow you specify in the `filter_values` parameter.
3. Use `include_md()` to reuse content fragments, adjusting them per workflow with `condition` and `replacements`.
     - add new content fragments to `common_text/`
     - make sure any **assets** in the included Markdown file are prefixed with `../../` to ensure they are correctly resolved relative to the final output location. This is important for images, files, and othrr assets to ensure they are correctly displayed in the final documentation. You may recieve a "WARNING" but this is expected and can be ignored if the assets are correctly displayed in the rendered documentation.

---

## For Developers

This section is for devs extending the macros or debugging their behavior.

### Macro Registration

Macros are defined in a `define_env(env)` function. This is the entry point that `mkdocs-macros-plugin` calls to register the Python logic.

### Macro Internals

#### `render_tsv_table(...)`

- Parses TSV files using `csv.DictReader`.
- Filters rows by matching values in a specific column.
- Expands wildcard patterns in column names.
- Sorts output by one or more columns.
- Outputs a standard Markdown table, with optional formatting for "input" rows.

#### `include_md(...)`

- Resolves a Markdown file path relative to `docs/`.
- Adjusts heading levels via `adjust_heading()`.
- Resolves relative links and image paths with `resolve_links()`.
- Detects and renders conditional blocks using `<!-- if:... -->` and `<!-- endif -->`.
- Recursively includes other Markdown files when `{{ include_md(...) }}` is found in the included file. Don't go too deep with this, as that just causes brain fog. All conditions and replacements are passed to the nested include as well.

### Supporting Functions

| Function                | Purpose |
|-------------------------|---------|
| `adjust_heading(...)`   | Normalizes Markdown heading levels with optional indentation. |
| `resolve_links(...)`    | Rewrites relative paths in `[]()` and `![]()` links to match the final structure. |
| `parse_macro_args(...)` | Parses nested macro arguments from Jinja-style include calls. |
