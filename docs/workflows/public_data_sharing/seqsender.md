# SeqSender

The SeqSender workflow submits samples from a Terra table to NCBI's BioSample, SRA and GenBank databases with [SeqSender](https://github.com/CDCgov/seqsender).

## Requirements

1. **A config file** with your NCBI account credentials and submission details.
2. **The required metadata columns in your Terra table.** A [validation-only run](#prepare-your-terra-table) lists them for you.

## Create Your Config File

Mirrored from SeqSender's Submission Wizard, this config builder runs in your browser and can be used to create a config file used in the SeqSender_PHB workflow. It will encrypt your NCBI FTP password and generate a corresponding key which will be required for any subsequent submissions. Note your key will change every time you create a new config file, even with the same password, unless you specify an **Existing key** beforehand.

--8<-- "docs/workflows/public_data_sharing/seqsender_config_builder.html"

Then, in your Terra workspace:

1. Upload the config file to your workspace's files and select it for the `seqsender_config` input.
2. Paste the key into the `seqsender_key` input.

## Prepare Your Terra Table

SeqSender reads each sample's metadata from your Terra table, so the table needs a column for every piece of metadata NCBI requires. Which columns those are depends on the databases you're submitting to and the organism you're working with.

1. **Run the workflow with `validation_only` set to `true`.** This checks every sample without submitting anything.
2. **Download the `validation_metadata_template` output.** It's a copy of your input Terra table, with empty columns added for each required piece of metadata. Every empty cell holds a description of what goes there, in angle brackets (`<...>`).
3. **Replace every `<...>` cell with your value.** Every cell with template text must be replaced, or it will fail validation. Open the file in Excel or Google Sheets, and keep `collection_date` as text in the form `YYYY-MM-DD` (spreadsheets like to reformat dates).
4. **Upload it to your Terra table** with **Import Data → Upload TSV** on the Data tab. Terra will update your existing table with the new metadata columns.
5. **Run with `validation_only` again** until every sample passes in `validation_summary`, then run without it to submit.


### Optional Columns

The template lists only the required columns. To send more metadata, add a column to your Terra table with the SeqSender column name, and the workflow passes it along. SeqSender's [metadata templates](https://github.com/CDCgov/seqsender/tree/v1.5.2/shiny/templates) list every column it accepts, with a description of each: `config.biosample.<your package>_template.csv` for BioSample, `config.sra.sra.schema_template.csv` for SRA, and `config.genbank.*` for GenBank.