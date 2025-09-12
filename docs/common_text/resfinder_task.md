
??? task "`ResFinder`: AMR Genotyping and XDR Shigella Prediction (optional)"
    To activate this task, set `call_resfinder` to `true`.

    **AMR Genotyping**

    The `ResFinder` task is an optional task that can be used in conjunction with AMRFinderPlus for detection and identification of AMR genes and resistance-associated mutations. This task runs the Centre for Genomic Epidemiology (CGE) ResFinder tool to identify acquired antimicrobial resistance. The default thresholds for calling AMR genes are 90% identity and 50% coverage of the reference genes (expressed as a fraction in workflow inputs: 0.9 and 0.5). These are the same thresholds utilized in BioNumerics for calling AMR genes.

    This task can also run the CGE PointFinder tool if the `call_pointfinder` variable is set with to `true`. The databases underlying the task are different to those used by AMRFinderPlus.

    ??? toggle "PointFinder-supported organisms"
        The following organisms are currently supported by _PointFinder_ for mutational-based predicted resistance. If GAMBIT ([see above](#taxonomic-assignment)) predicted the sample to be an organism not on this list, PointFinder will be skipped.
       
        - Campylobacter coli & C. jejuni
        - Enterococcus faecalis
        - Enterococcus faecium
        - Escherichia coli & Shigella spp.
        - Helicobacter pylori
        - Neisseria gonorrhoeae
        - Klebsiella
        - Mycobacterium tuberculosis
        - Salmonella spp.
        - Staphylococcus aureus

    **XDR Shigella prediction**

    The `ResFinder` Task also has the ability to predict whether or not a sample meets the CDC's definition for extensively drug-resistant (XDR) Shigella. 

    > _CDC defines XDR Shigella bacteria as strains that are resistant to all commonly recommended empiric and alternative antibiotics — azithromycin, ciprofloxacin, ceftriaxone, trimethoprim-sulfamethoxazole (TMP-SMX), and ampicillin_. [See also the _Increase in Extensively Drug-Resistance Shigellosis in the United State_ CDC Health Network Alert](https://emergency.cdc.gov/han/2023/han00486.asp) where this definition can be found.
    
    ??? toggle "Criteria for XDR Shigella Prediction"
        A sample is required to meet **all 7 criteria** in order to be predicted as `XDR Shigella` 

        1. The GAMBIT task in the workflow must identify the sample as `Shigella` OR the user must input the word `Shigella` somewhere within the input String variable called `expected_taxon`. This requirement serves as the identification of a sample to be of the Shigella genus.
        2. Resfinder or PointFinder predicted resistance to **Ampicillin**
        3. Resfinder or PointFinder predicted resistance to **Azithromycin**
        4. Resfinder or PointFinder predicted resistance to **Ciprofloxacin**
        5. Resfinder or PointFinder predicted resistance to **Ceftriazone**
        6. Resfinder or PointFinder predicted resistance to **Trimethoprim**
        7. Resfinder or PointFinder predicted resistance to **Sulfamethoxazole**

    There are 3 potential outputs for the **`resfinder_predicted_xdr_shigella`** output string**:**

    - **`Not Shigella based on gambit_predicted_taxon or user input`**
    - **`Not XDR Shigella`** for samples identified as Shigella by GAMBIT or user input BUT does ResFinder did not predict resistance to **all 6 drugs in XDR definition**
    - **`XDR Shigella`** meaning the sample was identified as Shigella and ResFinder/PointFinder did predict resistance to ceftriazone, azithromycin, ciprofloxacin, trimethoprim, sulfamethoxazole, and ampicillin.
    
    !!! techdetails "ResFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_resfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_resfinder.wdl) |
        | Software Source Code | [ResFinder Tool on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder/src/master/)<br>[ResFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)<br>[PointFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/) |
        | Software Documentation | [ResFinder on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder/src/master/README.md) |
        | Original Publication(s) | _ResFinder tool_: [ResFinder 4.0 for predictions of phenotypes from genotypes](https://academic.oup.com/jac/article/75/12/3491/5890997)<br>_ResFinder database_: [Identification of acquired antimicrobial resistance genes](https://doi.org/10.1093/jac/dks261) |
