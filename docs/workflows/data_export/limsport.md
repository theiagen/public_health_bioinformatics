# LIMSport

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filters={"Name": "[**LIMSport**](../workflows/data_export/limsport.md)"}, columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level", "Dockstore"]) }}

## LIMSport_PHB

[LIMSport](https://github.com/theiagen/limsport) generates a LIMS-importable file, in addition to a QC fail report and a QC summary file. It can be used to replace hand-curation of Terra results before LIMS import.

In order, LIMSport_PHB:

1. downloads a Terra table,
2. extracts the specified samples in `sample_names`,
3. selects all columns that are indicated in the configuration file,
4. applies QC to each row according to that config, then
5. reorders and renames the columns for the output report

This workflow runs at the set-level and produces a single report for all listed samples.

### Configuration File

The configuration file decides everything about the output file. Please see the [LIMSport config documentation](https://theiagen.github.io/limsport/latest/#the-config-file) for more details.

Briefly, the configuration file is a YAML file with two primary keys, `sample_qc` and `set_qc`, with an optional third key `skip_sample_qc`.

???+ dna "`sample_qc`"
    `sample_qc` lists every column to include in the LIMS report and in the order they appear. Each column can be renamed, be used only for QC, or be written to the output report after passing QC or without QC altogether. The sample ID column (the first column of the input table) is always written, even if it isn't listed.

    QC conditions can be either fixed thresholds (e.g., the value of column X must be greater than the threshold Y) or conditional thresholds (e.g., a different genome-length range per predicted taxon). Rows that do not match any conditional value fail unless a default is set.

???+ dna_blue "`set_qc`"
    `set_qc` checks specific sample(s) for specified metrics, and if any of its checks fail, the entire run will fail QC, resulting in an empty LIMS report. This is useful for gauging a run's quality based on control samples, like negative controls or positive controls not performing as expected.

    Please note that all samples indicated in `set_qc` rules must be included in `sample_names`. If not, the workflow will fail.

???+ dna "`skip_sample_qc`"
    `skip_sample_qc` can be used to indicate which particular samples should skip per-row QC. This is most useful for negative controls, where the samples are still checked by `set_qc` but regular QC conditions do not apply. These samples are marked as "SKIPPED" in the `sample_qc_status` column in the QC summary output, and are included in the final LIMS report.

#### Config Builder

Use the following config builder to build your configuration file from scratch, or drag and drop an existing one to make edits.

<!-- The builder is served from LIMSport's own GitHub Pages site, so PHB never copies its code. Keep the version in this URL matching the limsport docker tag in tasks/utilities/data_handling/task_limsport.wdl. `latest` follows the newest LIMSport release; swap to a pinned version (e.g. v0.1.2) once that docs version is deployed.
allow="clipboard-write" is required for the builder's "Copy to clipboard" button: browsers block clipboard access from inside an iframe served by another website unless the parent page grants it explicitly. data-sync-theme makes docs/theme/javascripts/iframe-theme.js send the site's light/dark scheme to the builder. -->
<iframe data-sync-theme src="https://theiagen.github.io/limsport/v1.0.0/config-builder/index.html" title="LIMSport Config Builder" width="100%" height="900" style="border: 1px solid var(--md-default-fg-color--lightest); border-radius: 4px;" allow="clipboard-write" loading="lazy"></iframe>

If the builder does not load above, [open it in a new tab](https://theiagen.github.io/limsport/latest/config-builder/index.html){:target="_blank"}.

### Inputs

!!! tip "QC fails are not included by default"
    The LIMS report will not include samples that fail QC unless **`include_qc_fails`** is set to true.

!!! warning "File parsing"
    File-parsing commands run only when **`allow_file_parsing`** is true; otherwise the workflow fails. **`max_file_parsing_failures`** fails the workflow once more than that many rows fail parsing (no limit by default).

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filters={"Workflow": "LIMSport"}, columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

### Workflow Tasks

{{ include_md("common_text/download_terra_table_task.md") }}

{{ include_md("common_text/limsport_task.md") }}

### Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filters={"Workflow": "LIMSport"}, columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

Please contact us if you need help writing a configuration file for a specific LIMS or your specific use-case!
