??? task "`Taxon Tables`: Copy outputs to new data tables based on taxonomic assignment (optional)"
    !!! tip ""
        This task is incompatible when running TheiaProk on the command-line as it is geared specifically for Terra. Do not activate this task if you are a command-line user.

    Activate this task by providing a value for the `taxon_tables` input variable. If provided, the user must also provide values to the `terra_project` and `terra_workspace` optional input variables.

    The `taxon_tables` module, if enabled, will copy sample data to a different data table based on the taxonomic assignment. For example, if an *E. coli* sample is analyzed, the module will copy the sample data to a new table for *E. coli* samples or add the sample data to an existing table.

    !!! tip ""
        To activate the `taxon_tables` module, provide a file indicating data table names to copy samples of each taxa to in the `taxon_tables` input variable.
        
        **Formatting the `taxon_tables` file**
        
        The `taxon_tables`  file must be uploaded a Google storage bucket that is accessible by Terra and should be in the format below. Briefly, the bacterial genera or species should be listed in the leftmost column with the name of the data table to copy samples of that taxon to in the rightmost column.
        
        | taxon | taxon_table |
        | --- | --- |
        | Listeria_monocytogenes | lmonocytogenes_specimen |
        | Salmonella | salmonella_specimen |
        | Escherichia | ecoli_specimen |
        | Shigella | shigella_specimen |
        | Streptococcus | strep_pneumo_specimen |
        | Legionella | legionella_specimen |
        | Klebsiella | klebsiella_specimen |
        | Mycobacterium | mycobacterium_specimen |
        | Acinetobacter | acinetobacter_specimen |
        | Pseudomonas | pseudomonas_specimen |
        | Staphylococcus | staphyloccus_specimen |
        | Neisseria | neisseria_specimen |
        
    !!! tip ""        
        There are no output columns for the taxon table task. The only output of the task is that additional data tables will appear for in the Terra workspace for samples matching a taxa in the `taxon_tables` file.

    !!! techdetails "`export_taxon_table` Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_export_taxon_table.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/utilities/data_export/task_export_taxon_table.wdl) |
