# Assembly Fetch

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Data Import](../../workflows_overview/workflows_type.md/#data-import) | [Any taxa](../../workflows_overview/workflows_kingdom.md/#any-taxa) | PHB v1.3.0 | Yes | Sample-level |

## Assembly_Fetch_PHB

The `Assembly_Fetch` workflow downloads assemblies from NCBI. This is particularly useful when you need to align reads against a reference genome, for example during a reference-based phylogenetics workflow. This workflow can be run in two ways:

1. You can provide an accession for the specific assembly that you want to download, and `Assembly_Fetch` will run only the NCBI genome download task to download this assembly,
2. You can provide an assembly, and `Assembly_Fetch` will first use the `ReferenceSeeker` task to first find the closest reference genome in RefSeq to your query assembly and then will run the NCBI genome download task to download that reference assembly.

!!! tip

    NOTE: If using Assembly_Fetch workflow version 1.3.0 or higher, the call-caching feature of Terra has been DISABLED to ensure that the workflow is run from the beginning and data is downloaded fresh. Call-caching will not be enabled, even if the user checks the box ✅ in the Terra workflow interface.

### Inputs

Assembly_Fetch requires the input samplename, and either the accession for a reference genome to download (ncbi_accession) or an assembly that can be used to query RefSeq for the closest reference genome to download (assembly_fasta).

This workflow runs on the sample level.

## Inputs

Assembly_Fetch requires the input samplename, and either the accession for a reference genome to download (ncbi_accession) or an assembly that can be used to query RefSeq for the closest reference genome to download (assembly_fasta).

This workflow runs on the sample level.
 
<div style="width: 100%;"> 
  <button id="resetOrderButton" style="padding: 4px 6px; font-size: 12px; background-color: #116eb7; color: white; border: none; border-radius: 5px; cursor: pointer;">
    Reset Table
  </button>

  <table id="inputTable" class="display" style="width: 100%; font-size: 14px;"> 
    <thead></thead>
    <tbody></tbody>
  </table>
</div>

<style>
  /* Ensure the table fills the container width */
  #inputTable {
    width: 100% !important;
    table-layout: auto;
  }

  /* Adjust font size for table headers and cells */
  #inputTable th, #inputTable td {
    font-size: 12px;
  }

  /* Style the table header */
  #inputTable th {
    background-color: #116eb7;
    color: white;
    font-weight: bold;
  }
</style>

<script>
  document.addEventListener('DOMContentLoaded', function() {
    function loadTable(enablePostbackSafe) {
      // Set default value if parameter is undefined
      if (typeof enablePostbackSafe === 'undefined') {
        enablePostbackSafe = true;
      }

      Papa.parse('/public_health_bioinformatics/data/assembly_fetch_input.csv', {
        download: true,
        header: true,
        complete: function(results) {
          var data = results.data;
          var headers = Object.keys(data[0]);  // Automatically fetch headers from CSV

          // Build table headers
          var headerHTML = '<tr>';
          headers.forEach(function(header) {
            headerHTML += '<th>' + header + '</th>';
          });
          headerHTML += '</tr>';
          document.querySelector('#inputTable thead').innerHTML = headerHTML;

          // Build table data from CSV
          var tableRows = [];
          data.forEach(function(row) {
            var rowArray = [];
            headers.forEach(function(header) {
              var cellContent = row[header] ? row[header] : '';
              if (header === 'Variable') {
                cellContent = '<strong>' + cellContent + '</strong>';
              }
              rowArray.push(cellContent);
            });
            tableRows.push(rowArray);
          });

          // Ensure DataTable is destroyed and reinitialized if it exists
          if ($.fn.DataTable.isDataTable('#inputTable')) {
            $('#inputTable').DataTable().clear().destroy();
          }

          // Initialize DataTable with scrolling options
          var table = $('#inputTable').DataTable({
            data: tableRows,
            columns: headers.map(function(header) {
              return { title: header };
            }),
            paging: false,
            searching: true,
            ordering: true,
            order: [],
            autoWidth: false,
            scrollY: "600px",  // Scrollable height
            scrollX: true,      // Enable horizontal scroll if needed
            scrollCollapse: true,  // Collapse the scroll area when fewer rows
            responsive: true,
            stateSave: true
          });

          // Adjust columns after initialization
          table.columns.adjust();

          // Initialize colResizable after DataTable
          $('#inputTable').colResizable({
            liveDrag: true,
            headerOnly: false,
            postbackSafe: enablePostbackSafe
          });
        }
      });
    }

    // Load the table when the DOM is ready
    loadTable();

    // Reset button to reload the table
    document.getElementById('resetOrderButton').addEventListener('click', function() {
      // Disable colResizable and flush its state
      $('#inputTable').colResizable({ disable: true, flush: true });

      // Clear DataTables state to reset column widths and other settings
      $('#inputTable').DataTable().state.clear();

      // Reload the table without saving column widths
      loadTable(false);
    });
  });
</script>



### Analysis Tasks

??? task "ReferenceSeeker (optional) Details"

    ##### ReferenceSeeker {#referenceseeker}

    `ReferenceSeeker` uses your draft assembly to identify closely related bacterial, viral, fungal, or plasmid genome assemblies in [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/).

    Databases for use with ReferenceSeeker are as follows, and can be used by pasting the gs uri in double quotation marks `" "` into the `referenceseeker_db` optional input:

    - archea:  `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-archaea-refseq-205.v20210406.tar.gz`
    - bacterial (**default**): `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-bacteria-refseq-205.v20210406.tar.gz`
    - fungi: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-fungi-refseq-205.v20210406.tar.gz`
    - plasmids: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-plasmids-refseq-205.v20210406.tar.gz`
    - viral: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-viral-refseq-205.v20210406.tar.gz`

    For ReferenceSeeker to identify a genome, it must meet user-specified thresholds for sequence coverage (`referenceseeker_conserved_dna_threshold`) and identity (`referenceseeker_ani_threshold`). The default values for these are set according to community standards (conserved DNA >= 69 % and ANI >= 95 %). A list of closely related genomes is provided in `referenceseeker_tsv`. The reference genome that ranks highest according to ANI and conserved DNA values is considered the closest match and will be downloaded, with information about this provided in the `assembly_fetch_referenceseeker_top_hit_ncbi_accession` output.

    !!! techdetails "ReferenceSeeker Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_referenceseeker.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_referenceseeker.wdl) |
        | Software version | 1.8.0 ("us-docker.pkg.dev/general-theiagen/biocontainers/referenceseeker:1.8.0--pyhdfd78af_0") |
        | Software Source Code | https://github.com/oschwengers/referenceseeker |
        | Software Documentation | https://github.com/oschwengers/referenceseeker |
        | Original Publication(s) | [ReferenceSeeker: rapid determination of appropriate reference genomes](https://joss.theoj.org/papers/10.21105/joss.01994) |

??? task "NCBI Datasets Details"

    ##### NCBI Datasets {#ncbi-datasets}

    The [`NCBI Datasets`](https://www.ncbi.nlm.nih.gov/datasets/) task downloads specified assemblies from NCBI using either the [virus](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/) or [genome](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/) (for all other genome types) package as appropriate.

    !!! techdetails "NCBI Datasets Technical Details"

        |  | Links |
        | --- | --- |
        | Task | [task_ncbi_datasets.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_import/task_ncbi_datasets.wdl) |
        | Software version | 14.13.2 (us-docker.pkg.dev/general-theiagen/staphb/ncbi-datasets:14.13.2) |
        | Software Source Code | https://github.com/ncbi/datasets |
        | Software Documentation | https://github.com/ncbi/datasets |
        | Original Publication(s) | Not known to be published |

### Outputs

| **Variable** | **Type** | **Description** |
|---|---|---|
| assembly_fetch_analysis_date | String | Date of assembly download |
| assembly_fetch_ncbi_datasets_assembly_data_report_json | File | JSON file containing report about assembly downloaded by Asembly_Fetch |
| assembly_fetch_ncbi_datasets_assembly_fasta | File | FASTA file downloaded by Assembly_Fetch |
| assembly_fetch_ncbi_datasets_docker | String | Docker file used for NCBI datasets |
| assembly_fetch_ncbi_datasets_gff | File | Assembly downloaded by Assembly_Fetch in GFF3 format |
| assembly_fetch_ncbi_datasets_gff3 | File | Assembly downloaded by Assembly_Fetch in GFF format |
| assembly_fetch_ncbi_datasets_version | String | NCBI datasets version used |
| assembly_fetch_referenceseeker_database | String | ReferenceSeeker database used |
| assembly_fetch_referenceseeker_docker | String | Docker file used for ReferenceSeeker |
| assembly_fetch_referenceseeker_top_hit_ncbi_accession | String | NCBI Accession for the top it identified by Assembly_Fetch |
| assembly_fetch_referenceseeker_tsv | File | TSV file of the top hits between the query genome and the Reference Seeker database |
| assembly_fetch_referenceseeker_version | String | ReferenceSeeker version used |
| assembly_fetch_version | String | The version of the repository the Assembly Fetch workflow is in |

## References

> **ReferenceSeeker:** Schwengers O, Hain T, Chakraborty T, Goesmann A. ReferenceSeeker: rapid determination of appropriate reference genomes. J Open Source Softw. 2020 Feb 4;5(46):1994.
<!-- -->
> **NCBI datasets: datasets:** NCBI Datasets is an experimental resource for finding and building datasets [Internet]. Github; [cited 2023 Apr 19]. Available from: <https://github.com/ncbi/datasets>
