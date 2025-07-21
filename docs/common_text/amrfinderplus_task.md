??? task "`AMRFinderPlus`: AMR Genotyping (default)"

    NCBI's [AMRFinderPlus](https://github.com/ncbi/amr/wiki) is the default antimicrobial resistance (AMR) detection tool used in TheiaProk. ResFinder may be used alternatively and if so, AMRFinderPlus is not run. 

    AMRFinderPlus identifies acquired antimicrobial resistance (AMR) genes, virulence genes, and stress genes.  Such AMR genes confer resistance to antibiotics, metals, biocides, heat, or acid. For some taxa (see [here](https://github.com/ncbi/amr/wiki/Running-AMRFinderPlus#--organism-option)), AMRFinderPlus will provide taxa-specific results including filtering out genes that are almost ubiquitous in the taxa (intrinsic genes) and identifying resistance-associated point mutations.  In TheiaProk, the taxon used by AMRFinderPlus is specified based on the `gambit_predicted_taxon` or a user-provided `expected_taxon`. AMRFinderPlus also has the ability to utilize a GFF and protein FASTA file which can be enabled via `amrfinder_use_gff` allowing for more accurate calls.

    You can check if a gene or point mutation is in the AMRFinderPlus database [here](https://www.ncbi.nlm.nih.gov/pathogens/refgene/#), find the sequences of reference genes [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047), and search the query Hidden Markov Models (HMMs) used by AMRFinderPlus to identify AMR genes and some stress and virulence proteins ([here](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)). The AMRFinderPlus database is updated frequently. You can ensure you are using the most up-to-date version by specifying the docker image as a workflow input. You might like to save this docker image as a workspace data element to make this easier.

    !!! techdetails "AMRFinderPlus Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_amrfinderplus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_amrfinderplus.wdl) |
        | Software Source Code | [amr on GitHub](https://github.com/ncbi/amr) |
        | Software Documentation | https://github.com/ncbi/amr/wiki |
        | Original Publication(s) | [AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/) |