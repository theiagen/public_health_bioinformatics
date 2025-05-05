# Terra_2_NCBI

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Public Data Sharing](../../workflows_overview/workflows_type.md/#public-data-sharing) | [Bacteria](../../workflows_overview/workflows_kingdom.md#bacteria), [Mycotics](../../workflows_overview/workflows_kingdom.md#mycotics) [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v3.0.0 | No | Set-level |

## Terra_2_NCBI_PHB

!!! warning "Do not resubmit!"
    **If the Terra_2_NCBI workflow fails, DO NOT resubmit.**

    Resubmission risks duplicate submissions and future failures.

    Contact Theiagen (`support@theiagen.com`) to determine the reason for failure, and **only move forward with Theiagen's guidance**.

!!! dna "Key Resources"
    - [Pathogen metadata formatter](../../assets/metadata_formatters/Terra_2_NCBI-PATHOGEN-metadata-2024-04-30.xlsx)
    - [Microbe metadata formatter](../../assets/metadata_formatters/Terra_2_NCBI-MICROBE-metadata-2022-07-11.xlsx)
    - [Virus metadata formatter](../../assets/metadata_formatters/Terra_2_NCBI-VIRUS-metadata-2022-09-09.xlsx)
    - [SARS-CoV-2 Wastewater metadata formatter](../../assets/metadata_formatters/Terra_2_NCBI-SC2WW-metadata-2025-01-10.xlsx)

The Terra_2_NCBI workflow is a programmatic data submission method to share metadata information with NCBI BioSample and paired-end Illumina reads with NCBI SRA directly from Terra without having to use the NCBI portal.

### Prerequisites

??? toggle "Before running the Terra_2_NCBI workflow"

    1. The user **must** have access to the NCBI FTP. To gain these credentials, we recommend emailing `**sra@ncbi.nlm.nih.gov**` a variation of the following example, including all the information:

        > Hello,
        >
        >We would like to automate submissions to the Submission Portal using XML metadata to accompany our cloud-hosted data files.  We would like to upload via FTP and need to create a submission group.
        >
        >Here is the relevant information:
        >
        >1. Suggested group abbreviation:
        >2. Full group name:
        >3. Institution and department:
        >4. Contact person (someone likely to remain at the location for an extended time):
        >5. Contact email:
        >6. Mailing address (including country and postcode):
        >
        >We will be using an existing submission pipeline that is known to work and would like to request that the production folder be activated. Thank you for your assistance!

    2. From NCBI, you will need to get in response:
        1. an FTP address (it will likely be ftp-private.ncbi.nih.gov)
        2. Username (typically the suggested group abbreviation)
        3. Password
        4. an acknowledgment that the production folder has been activated.

        Please confirm that the production folder has been activated, or else the submission pipeline will either fail or only run test submissions and not actually submit to NCBI.

    3. Before you can run the workflow for the first time, we also recommend scheduling a meeting with Theiagen to get additional things set up, including
        - adding a correctly-formatted configuration file to your workspace data elements that includes your FTP username and password, laboratory details, and other important information.
        - ensuring your proxy account has been given permission to write to the google bucket where SRA reads are temporarily stored before being transferred to NCBI.

        ??? toggle "What is the configuration file used for?"
            The configuration file tells the workflow your username and password so you can access the FTP. It also provides important information about who should be contacted regarding the submission. We recommend contacting a member of Theiagen for help in the creation of this configuration file to ensure that everything is formatted correctly.

### Collating BioSample Metadata

In order to create BioSamples, you need to choose the correct BioSample package and have the appropriate metadata included in your data table.

Currently, Terra_2_NCBI only supports _Pathogen_, _Virus_, _Microbe_, and _SARS-CoV-2 Wastewater Surveillance_ BioSample packages. **Most organisms should be submitted using the Pathogen package** unless you have been specifically directed otherwise (either through CDC communications or another reliable source). Definitions of packages supported by Terra_2_NCBI are listed below with more requirements provided via the links:

- [Pathogen.cl](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Pathogen.cl.1.0/) - any clinical or host-associated pathogen
- [Pathogen.env](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Pathogen.env.1.0/) - environmental, food or other pathogen *(no metadata formatter available at this time)*
- [Microbe](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Microbe.1.0/) - bacteria or other unicellular microbes that do not fit under the MIxS, Pathogen, or Virus packages.
- [Virus](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/Virus.1.0/) -  viruses **not** directly associated with disease
    - Viral pathogens should be submitted using the Pathogen: Clinical or host-associated pathogen package.
- [SARS-CoV-2.wwsurv](https://www.ncbi.nlm.nih.gov/biosample/docs/packages/SARS-CoV-2.wwsurv.1.0/) - SARS-CoV-2 wastewater surveillance samples

### Metadata Formatters

For each package, we have created a metadata template spreadsheet to help you organize your metadata:

Please note that the pathogen metadata formatter is for the _clinical_ pathogen package, not the environmental pathogen.

- [Terra_2_NCBI-PATHOGEN-metadata-2024-04-30.xlsx](../../assets/metadata_formatters/Terra_2_NCBI-PATHOGEN-metadata-2024-04-30.xlsx)
- [Terra_2_NCBI-MICROBE-metadata-2022-07-11.xlsx](../../assets/metadata_formatters/Terra_2_NCBI-MICROBE-metadata-2022-07-11.xlsx)
- [Terra_2_NCBI-VIRUS-metadata-2022-09-09.xlsx](../../assets/metadata_formatters/Terra_2_NCBI-VIRUS-metadata-2022-09-09.xlsx)
- [Terra_2_NCBI-SC2WW-metadata-2025-01-10.xlsx](../../assets/metadata_formatters/Terra_2_NCBI-SC2WW-metadata-2025-01-10.xlsx)

We are constantly working on improving these spreadsheets and they will be updated in due course.

### Running the Workflow

We recommend running a test submission before your first production submission to ensure that all data has been formatted correctly. Please contact Theiagen (<support@theiagen.com>) to get this set up.

In the test submission, any real BioProject accession numbers you provide will not be recognized. You will have to make a "fake" or "test" BioProject. This cannot be done through the NCBI portal. Theiagen can provide assistance in creating this as it requires manual command-line work on the NCBI FTP using the account they provided for you.

??? toggle "**What's the difference between a test submission and a production submission?**"

    A production submission means that your submission using Terra_2_NCBI will be submitted to NCBI as if you were using the online portal. That means that anything you submit on production will be given to the ****real** **NCBI servers and appear and become searchable on the NCBI website.
    
    A test submission gives your data to a completely detached **replica** of the production server. This means that any data you submit as a test will behave exactly **like a real submission, but since it's detached, nothing **will appear on the NCBI website, and anything returned from the workflow (such as BioSample accession numbers) will be fake. If you search for these test BioSample accession numbers on the NCBI website, either (a) nothing will appear, or (b) it will link to a random sample. 
    
    If you want your data to be on NCBI, you must run a production submission. Initially, NCBI locks the production folder so that the user doesn't accidentally submit test data to the main database. You must have requested activation of the production folder prior to your first production submission.

### Inputs

This workflow runs on set-level data tables.

!!! info "Production Submissions"
    Please note that an optional Boolean variable, `submit_to_production`, is **required** for a production submission.

???+ tip "Using Customized Column Names in Terra Tables"

    In some cases, users may have data tables in Terra with column names that differ from the default expected by the workflow. The `Terra_2_NCBI` workflow allows users to supply a **custom column mapping file**, enabling them to specify how their columns map to the required workflow variables.

    To use a custom column mapping file:

    1. Create a tab-delimited `.tsv` file with the following structure:
   
        A header including "Custom" and "Required" should be included in the first row.
        The "Custom" column should contain the actual column names in your Terra table (e.g., 'collection-date'), and the "Required" column should contain the column names expected by the workflow (e.g., 'collection_date').

        Example Mapping File:
        ```plaintext
        Custom  Required
        Collection-Date collection_date
        geo_location    geo_loc_name
        bioproject_column   bioproject
        sample_id_column    sample_names
        ```

    2. Upload the file to your Terra workspace and reference it in the `column_mapping_file` parameter when running the workflow using Google Cloud Storage paths.

    Ensure the mapping file includes all columns with custom names. Columns that match the default workflow names do not need to be included. Missing mappings for renamed columns may result in errors during execution if the column is required, and will not be found if the column is optional.

    The workflow will automatically map the specified column names from your Terra table to the required workflow variables using the 'custom_mapping_file'.

    To find a list of the expected required and optional column names, [please refer to the code blocks that can be found here](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/submission/task_submission.wdl#L65). The required and optional metadata fields are organized by the BioSample type.

    Below, you can find the required metadata fields for the currently supported BioSample types:

    ??? toggle "Microbe Required Metadata"
        - submission_id
        - organism
        - collection_date
        - geo_loc_name
        - sample_type

    ??? toggle "Wastewater Required Metadata"
        - submission_id
        - organism
        - collection_date
        - geo_loc_name
        - isolation_source
        - ww_population
        - ww_sample_duration
        - ww_sample_matrix
        - ww_sample_type
        - ww_surv_target_1
        - ww_surv_target_1_known_present

    ??? toggle "Pathogen.cl Required Metadata"
        - submission_id
        - organism
        - collected_by
        - collection_date
        - geo_loc_name
        - host
        - host_disease
        - isolation_source
        - lat_lon
    
    ??? toggle "Pathogen.env Required Metadata"
        - submission_id
        - organism
        - collected_by
        - collection_date
        - geo_loc_name
        - isolation_source
        - lat_lon

    ??? toggle "Virus Required Metadata"
        - submission_id
        - organism
        - isolate
        - collection_date
        - geo_loc_name
        - isolation_source

    ---
    For further assistance in setting up a custom column mapping file, please contact Theiagen at [support@theiagen.com](mailto:support@theiagen.com).

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/input_tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Terra_2_NCBI", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

??? task "Workflow Tasks"

    ##### Workflow Tasks {% raw %} {#workflow-tasks} {% endraw %}

    The workflow will perform the following tasks, each highlighted as `code`

    1. `prune_table`formats all incoming metadata for submission.
    2. If you are submitting BioSamples:
        1. `biosample_submit_tsv_ftp_upload` will 
            1. format the BioSample table into XML format
            2. submit BioSamples to NCBI
            3. return all NCBI communications in XML format, and
            4. parse those communications for any and all BioSample accessions.
        2. `add_biosample_accessions` will
            1. add the BioSample accessions to SRA metadata
            2. upload the BioSample accessions to the origin Terra table

            If BioSample accessions fail to be generated, this task ends the workflow and users should contact Theiagen for further support. Otherwise, the workflow will continue and outputs are returned to the Terra data table. 

    3. If BioSample accessions were generated or if BioSample submission was skipped
        1. `sra_tsv_to_xml` converts the SRA metadata (including any generated or pre-provided BioSample accessions) into XML format.
        2. `ncbi_sftp_upload` 
            1. uploads the SRA metadata to NCBI 
            2. returns any XML communications from NCBI.

#### Workflow Success

If the workflow ends successfully, it returns the outputs to the Terra data table and the XML communications from NCBI will say that submission is underway. The workflow does not declare successful sample submission since SRA sometimes takes a while to do this. If the submission was successful, the point of contact for the submission will receive the SRA accessions via email from NCBI.

If the workflow ends unsuccessfully, no outputs will be shown on Terra and the `biosample_status` output variable will indicate that the BioSample submission failed.

### Outputs

The output files contain information mostly for debugging purposes. Additionally, if your submission is successful, the point of contact for the submission should also receive an email from NCBI notifying them of their submission success.

/// html | div[class="searchable-table"]

{{ input_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Terra_2_NCBI", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

???+ toggle "An example excluded_samples.tsv file"

    ##### An example excluded_samples.tsv file {% raw %} {#example-excluded-samples} {% endraw %}

    Due to the nature of tsv files, it may be easier to download and open this file in Excel. 
    
    [example_excluded_samples.tsv](../../assets/files/example_excluded_samples.tsv)
    
    ```
    Samples excluded for quality thresholds:
    sample_name message 
    sample2 VADR skipped due to poor assembly
    sample3 VADR number alerts too high: 3 greater than limit of 0
    sample4 Number of Ns was too high: 10000 greater than limit of 5000
    
    Samples excluded for missing required metadata (will have empty values in indicated columns):
    tablename_id    organism    country library_layout
    sample5         paired
    sample6 SARS-CoV-2  USA
    ```
    
    This example informs the user that samples 2-4 were excluded for quality reasons (the exact reason is listed in the `message` column), and that samples 5 and 6 were excluded because they were missing required metadata fields (sample5 was missing the `organism` and `country` fields, and sample6 was missing the `library_layout` field).

### Limitations

- The maximum number of samples that can be submitted at once appears to be 300. We recommend submitting less than 300 samples at a time to avoid errors due to large submission sizes.
- A workflow on returning SRA accessions using the generated BioSample accessions is in progress.

### Acknowledgments

This workflow would not have been possible without the invaluable contributions of Dr. Danny Park.
