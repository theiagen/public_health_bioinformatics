# Terra_2_ENA

## Quick Facts

{{ render_tsv_table("docs/assets/tables/all_workflows.tsv", sort_by="Name", filter_column="Name", filter_values="[**Terra_2_ENA**](../workflows/public_data_sharing/terra_2_ena.md)", columns=["Workflow Type", "Applicable Kingdom", "Last Known Changes", "Command-line Compatibility","Workflow Level"]) }}

## Terra_2_ENA_PHB

This workflow utilizes the [ENA Webin-CLI Bulk Submission Tool](https://github.com/enasequence/ena-bulk-webincli) to bulk submit read data to the ENA.

## ENA Submissions

Before you can submit data to ENA you must [register](https://www.ebi.ac.uk/ena/submit/webin/login) a Webin submission account. ENA allows submissions via the Webin-CLI program which validates your submissions entirely before you complete them. Submissions made through Webin are represented using a number of different [metadata objects](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/metadata.html). Submissions to ENA result in [accession numbers](https://ena-docs.readthedocs.io/en/latest/submit/general-guide/accessions.html) and these accessions can be used to identify each unique part of your submission. See the [ENA submission documentation](https://ena-docs.readthedocs.io/en/latest/submit/general-guide.html) for more information.

## Pre-requisites

!!! warning "Before running the `Terra_2_ENA` workflow, make sure you have registered a study and are using the correct study accession number"

- To submit data into ENA you must first register a study to contain and manage it. Studies (also referred to as projects) can be registered through the Webin Portal. [Log in](https://www.ebi.ac.uk/ena/submit/webin/login) with your Webin credentials and select the ‘Register Study’ button to bring up the interface. Once registration is complete, you will be assigned accession numbers. You may return to the dashboard and select the ‘Studies Report’ button to review registered studies.

- Additionally, before submitting most types of data to ENA, samples must be [registered](https://ena-docs.readthedocs.io/en/latest/submit/samples.html#). To register samples, ensure that your Terra data table includes all the samples you intend to submit, along with their raw read data (`FASTQ`, `BAM`, or `CRAM` format) and associated metadata. To meet ENA’s requirements, each sample must include a minimum set of metadata. See below for the mandatory and recommended metadata fields, as well as the default column names used to identify them in your Terra data table.

## What needs to be included in your Terra data table?

### Read Data Fields

<div class="grid cards" markdown>

-   ???+ toggle "**Mandatory Fields**"

        !!! warning "These columns are required for submission and must be included in the Terra data table. The column names must appear exactly as shown and cannot be substituted or modified using column mappings."

        | Terra Column Name | Description |
        |---|---|
        | `read1`/`read2`/`bam_file`/`cram_file` | The path to two paired end `FASTQ` files, `BAM` file, or `CRAM` file containing sequencing data. |
        | `experiment_name` | Unique name of the experiment. |
        | `sequencing_platform` | The platform used to generate the sequence data. [See permitted values](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#platform). |
        | `sequencing_instrument` | The instrument used to generate the sequence data. [See permitted values](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#instrument). |
        | `library_source` | The source of the library. [See permitted values](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#source). |
        | `library_selection` | The method used to select the library. [See permitted values](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#selection). |
        | `library_strategy` | The strategy used to generate the library. [See permitted values](https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#strategy). |

</div>

<div class="grid cards" markdown>

-   ???+ toggle "**Optional Fields**"

        | Terra Column Name | Description |
        |---|---|
        | `insert_size` | The insert size for paired reads. |
        | `library_description` | Free text library description. |

</div>

### Sample Metadata Fields
??? tip "Using Customized Column Names in Terra Tables"

    In some cases, users may have data tables in Terra with column names that differ from the field names expected by ENA. The `Terra_2_ENA` workflow allows users to supply a **custom column mapping file**, enabling them to specify how their columns map to the required/mandatory field names.

    To use a custom column mapping file:

    1. Create a tab-delimited `.tsv` file with the following structure:

        A header including `terra_column` and `ena_column` should be included in the first row.
        The `terra_column` column should contain the actual column names in your Terra table (e.g., 'my_fav_collection_date'), and the `ena_column` column should contain the column names expected by ENA (e.g., `collection date`). More information about the mandatory and recommended metadata fields and associated column names are described in the tables below.

        Example Mapping File:
        ```plaintext
        terra_column            ena_column
        my_sample_title         title
        my_fav_collection_date  collection_date
        my_geo_loc_name         geo_loc_name
        ```

    2. Upload the file to your Terra workspace and reference it in the `column_mappings` input parameter when running the workflow using Google Cloud Storage paths.

    Ensure the mapping file includes all columns with custom names. Columns that match the default workflow names do not need to be included. Missing mappings for renamed columns may result in errors during execution if the column is required, and will not be found if the column is optional. The workflow will automatically map the specified column names from your Terra table to the required ENA field names as long as the mapping file is provided correctly.

=== "Bacterial Metadata"

    <div class="grid cards" markdown>

    -   ???+ toggle "**Mandatory Fields**"

            !!! warning "These fields are required for submission and must be included in the Terra data table or supplied as an input parameter"

                If you cannot provide a value for a mandatory field within, set the `allow-missing` input parameter to `true` or alternatively, use one of the [INDSC accepted terms](https://ena-docs.readthedocs.io/en/latest/submit/samples/missing-values.html) for missing value reporting.

            | <div style="width: 160px;">Terra Column Name</div> | <div style="width: 140px;">ENA Field Name</div> | Description |
            |---|---|---|
            | `title` | sample_title  | Title of the sample. |
            | `taxon_id` and/or `organism` | tax_id and/or scientific name  | Taxonomic identifier (NCBI taxon ID) or scientific name of the organism from which the sample was obtained. |
            | `collection_date` | collection date | The date the sample was collected with the intention of sequencing, either as an instance (single point in time) or interval. In case no exact time is available, the date/time can be right truncated i.e. all of these are valid ISO8601 compliant times: 2008-01-23T19:23:10+00:00; 2008-01-23T19:23:10; 2008-01-23; 2008-01; 2008. |
            | `geo_loc_name` | geographic location (country and/or sea) | The geographical origin of where the sample was collected from, with the intention of sequencing, as defined by the country or sea name. Country or sea names should be chosen from the INSDC country [list](http://insdc.org/country.html). |
            | `host_health_state` | host health state | Health status of the host at the time of sample collection. Must be one of the following: `diseased`, `healthy`, `missing: control sample`, `missing: data agreement established pre-2023`, `missing: endangered species`, `missing: human-identifiable`, `missing: lab stock`, `missing: sample group`, `missing: synthetic construct`, `missing: third party data`, `not applicable`, `not collected`, `not provided`, `restricted access`. |
            | `host_scientific_name` | host scientific name | Scientific name of the natural (as opposed to laboratory) host to the organism from which sample was obtained. |
            | `isolation_source` | isolation_source | Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived. |
            | `isolate` | isolate | Individual isolate from which the sample was obtained. |

    </div>

    <div class="grid cards" markdown>

    -   ???+ toggle "**Optional Fields**"

            | <div style="width: 190px;">Terra Column Name</div> | <div style="width: 190px;">ENA Field Name</div> | Description |
            |---|---|---|
            | `library_description` | sample_description | Description of the sample. |
            | `lat_lon` | lat_lon | Geographical coordinates of the location where the specimen was collected. |
            | `serovar` | serovar | Serological variety of a species (usually a prokaryote) characterized by its antigenic properties. |
            | `strain` | strain | Name of the strain from which the sample was obtained. |

    </div>

    Reference: [ENA prokaryotic pathogen minimal sample checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000028)

=== "Viral Metadata"

    <div class="grid cards" markdown>

    -   ???+ toggle "**Mandatory Fields**"

            !!! warning "These fields are required for submission and must be included in the Terra data table or supplied as an input parameter"

                If you cannot provide a value for a mandatory field within, set the `allow-missing` input parameter to `true` or alternatively, use one of the [INDSC accepted terms](https://ena-docs.readthedocs.io/en/latest/submit/samples/missing-values.html) for missing value reporting.

            | <div style="width: 170px;">Terra Column Name</div> | <div style="width: 140px;">ENA Field Name</div> | Description |
            |---|---|---|
            | `title` | sample_title  | Title of the sample. |
            | `taxon_id` and/or `organism` | tax_id and/or scientific name  | Taxonomic identifier (NCBI taxon ID) or scientific name of the organism from which the sample was obtained. |
            | `collection_date` | collection date | The date the sample was collected with the intention of sequencing, either as an instance (single point in time) or interval. In case no exact time is available, the date/time can be right truncated i.e. all of these are valid ISO8601 compliant times: 2008-01-23T19:23:10+00:00; 2008-01-23T19:23:10; 2008-01-23; 2008-01; 2008. |
            | `collecting_institution` | collecting institution | Name of the institution to which the person collecting the specimen belongs. Format: Institute Name, Institute Address |
            | `collector_name` | collector name | Name of the person who collected the specimen. Example: John Smith |
            | `geo_loc_name` | geographic location (country and/or sea) | The geographical origin of where the sample was collected from, with the intention of sequencing, as defined by the country or sea name. Country or sea names should be chosen from the INSDC country [list](http://insdc.org/country.html). |
            | `host_common_name` | host common name | Common name of the host, e.g. human |
            | `host_health_state` | host health state | Health status of the host at the time of sample collection. Must be one of the following: `diseased`, `healthy`, `missing: control sample`, `missing: data agreement established pre-2023`, `missing: endangered species`, `missing: human-identifiable`, `missing: lab stock`, `missing: sample group`, `missing: synthetic construct`, `missing: third party data`, `not applicable`, `not collected`, `not provided`, `restricted access`. |
            | `host_scientific_name` | host scientific name | Scientific name of the natural (as opposed to laboratory) host to the organism from which sample was obtained. |
            | `host_sex` | host sex | Gender or sex of the host. Must be one of the following: `female`, `male`, `hermaphrodite`, `neuter`, `not applicable`, `not collected`, `not provided`, `other`, `missing: control sample`, `missing: data agreement established pre-2023`, `missing: endangered species`, `missing: human-identifiable`, `missing: lab stock`, `missing: sample group`, `missing: synthetic construct`, `missing: third party data`. |
            | `host_subject_id` | host subject id | A unique identifier by which each subject can be referred to, de-identified, e.g. #131 |
            | `isolation_source` | isolation_source | Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived. |
            | `isolate` | isolate | Individual isolate from which the sample was obtained. |

    </div>

    <div class="grid cards" markdown>

    -   ???+ toggle "**Optional Fields**"

            | <div style="width: 190px;">Terra Column Name</div> | <div style="width: 190px;">ENA Field Name</div> | Description |
            |---|---|---|
            | `latitude` | geographic location (latitude) | The geographical origin of the sample as defined by latitude. The values should be reported in decimal degrees and in WGS84 system |
            | `longitude` | geographic location (longitude) | The geographical origin of the sample as defined by longitude. The values should be reported in decimal degrees and in WGS84 system |
            | `region_locality` | geographic location (region and locality) | The geographical origin of the sample as defined by the specific region name followed by the locality name. |
            | `host_disease_outcome` | host disease outcome | Disease outcome in the host. |
            | `host_age` | host age | Age of host at the time of sampling; relevant scale depends on species and study, e.g. could be seconds for amoebae or centuries for trees |
            | `host_behaviour` | host behaviour | Natural behaviour of the host. |
            | `host_habitat` | host habitat | Natural habitat of the avian or mammalian host. |
            | `isolation_source_host` | isolation source host-associated | Name of host tissue or organ sampled for analysis. Example: tracheal tissue |
            | `isolation_source_non_host` | isolation source non-host-associated | Describes the physical, environmental and/or local geographical source of the biological sample from which the sample was derived. Example: soil |
            | `receipt_date` | receipt date | Date on which the sample was received. Format:YYYY-MM-DD. Please provide the highest precision possible. If the sample was received by the institution and not collected, the 'receipt date' must be provided instead. |
            | `sample_capture_status` | sample capture status | Reason for the sample collection. |
            | `library_description` | sample_description | Description of the sample. |
            | `serotype` | serotype | Serological variety of a species characterised by its antigenic properties. For Influenza, HA subtype should be the letter H followed by a number between 1-16 unless novel subtype is identified and the NA subtype should be the letter N followed by a number between 1-9 unless novel subtype is identified. |
            | `virus_identifier` | virus identifier | Unique laboratory identifier assigned to the virus by the investigator. Strain name is not sufficient since it might not be unique due to various passsages of the same virus. Format: up to 50 alphanumeric characters |

    </div>

    Reference: [ENA viral minimal sample checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000033)


## Workflow Inputs

It's important to note that the `Terra_2_ENA` workflow is designed to run on set-level data tables. This means that the workflow will process all samples within a set together, rather than handling each sample individually. The `samples` input variable expects an array of sample IDs, corresponding to a set table. In most cases, set tables are generated automatically when running a workflow. However, if you need to create one manually, refer to this guide on [how to create a set table](https://support.terra.bio/hc/en-us/articles/6660506445339-How-to-make-a-set-table).

!!! warning "The `submit_to_production` input parameter is set to `false` by default. This means that the workflow will not submit data to the production ENA server unless you explicitly set it to `true`. This is useful for testing purposes, allowing you to validate your data without making actual submissions."

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_inputs.tsv", input_table=True, filter_column="Workflow", filter_values="Terra_2_ENA", columns=["Terra Task Name", "Variable", "Type", "Description", "Default Value", "Terra Status"], sort_by=[("Terra Status", True), "Terra Task Name", "Variable"]) }}

///

## Workflow Outputs

/// html | div[class="searchable-table"]

{{ render_tsv_table("docs/assets/tables/all_outputs.tsv", input_table=False, filter_column="Workflow", filter_values="Terra_2_ENA", columns=["Variable", "Type", "Description"], sort_by=["Variable"]) }}

///

## References
- [ENA prokaryotic pathogen minimal sample checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000028)
- [ENA viral minimal sample checklist](https://www.ebi.ac.uk/ena/browser/view/ERC000033)
- [ENA Data Submission Documentation](https://ena-docs.readthedocs.io/en/latest/index.html)