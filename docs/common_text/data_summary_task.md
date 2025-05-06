
??? task "Data summary (optional)"
    ##### Data Summary (optional)

    If you fill out the `data_summary_*` and `sample_names` optional variables, you can use the optional `summarize_data` task. The task takes a comma-separated list of column names from the Terra data table, which should each contain a list of comma-separated items. For example, `"amrfinderplus_virulence_genes,amrfinderplus_stress_genes"` (with quotes, comma separated, no spaces) for these output columns from running TheiaProk. The task checks whether those comma-separated items are present in each row of the data table (sample), then creates a CSV file of these results. The CSV file indicates presence (TRUE) or absence (empty) for each item. By default, the task does not add a Phandango coloring tag to group items from the same column, but you can turn this on by setting `phandango_coloring` to `true`.

    ??? toggle "**Example output CSV**"

        ```text linenums="1"
        Sample_Name,aph(3')-IIa,blaCTX-M-65,blaOXA-193,tet(O)
        sample1,TRUE,,TRUE,TRUE
        sample2,,,FALSE,TRUE
        sample3,,,FALSE,
        ```

    ??? toggle "**Example use of Phandango coloring**"

        Data summary produced using the `phandango_coloring` option, visualized alongside Newick tree at <http://jameshadfield.github.io/phandango/#/main>

        !!! caption "Example phandango_coloring output"
            ![Phandango coloring example](../../assets/figures/example_phandango_coloring.png)

    !!! techdetails "Data summary technical details"

        |  | Links |
        | --- | --- |
        | Task | [task_summarize_data.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_handling/task_summarize_data.wdl) |
