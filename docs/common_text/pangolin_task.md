??? task "`pangolin`"

    Pangolin (Phylogenetic Assignment of Named Global Outbreak Lineages) was developed to implement a dynamic nomenclature for designating SARS-CoV-2 lineage assignments and is used by researchers and public health agencies worldwide to track the spread and transmission of SARS-CoV-2.
    
    Pangolin aligns input sequences against an early SARS-CoV-2 reference and generates a unique hash for the alignment. The hash is checked against a designation cache to see if it matches any previously identified lineages, and is checked via [scorpio](https://github.com/cov-lineages/scorpio) (_Serious Constellations of Reoccuring Phylogenetically-Independent Origin_) to determine if the hash matches any variant of concern (VOC) [_constellations_](https://github.com/cov-lineages/constellations/tree/main), which are groups of functionally meaningful mutations that can independently evolve. Following a QC check, an inference pipeline is run: either [pangoLEARN](https://cov-lineages.org/resources/pangolin/pangolearn.html) or [UShER](https://github.com/yatisht/usher) (which is the default inference model). The final lineage report is then generated.

    !!! techdetails "Pangolin Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_pangolin.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/betacoronavirus/task_pangolin.wdl) |
        | Software Source Code | [Pangolin on GitHub](https://github.com/cov-lineages/pangolin) |
        | Software Documentation | [Pangolin on cov-lineages.org](https://cov-lineages.org/resources/pangolin.html) |
        | Original Publication(s) | [A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology](https://doi.org/10.1038/s41564-020-0770-5) |
