
??? task "`ResFinder`: AMR Genotyping and XDR Shigella Prediction (optional)"
    To activate this task, set `call_resfinder` to `true`.

    **AMR Genotyping**

    The `ResFinder` task is an optional task that can be used in conjunction with AMRFinderPlus for detection and identification of AMR genes and resistance-associated mutations. This task runs the Centre for Genomic Epidemiology (CGE) ResFinder tool to identify acquired antimicrobial resistance, and uses different databases for AMR detection than AMRFinderPlus. In order to call an AMR gene, there must be at least 90% (0.9) sequence identity and 50% (0.5) coverage of the reference genes. These are the same thresholds used in BioNumerics for AMR detection. 

    Set `call_pointfinder` to `true` to also run the CGE PointFinder tool for detection of chromosomal mutations associated with antimicrobial resistance. 

    ??? toggle "PointFinder-supported organisms"
        The following organisms are currently supported by _PointFinder_ for mutational-based predicted resistance. If GAMBIT ([see above](#taxonomic-assignment)) predicted the sample to be an organism not on this list, PointFinder will be skipped.
       
        - Campylobacter coli & C. jejuni
        - Campylobacter spp. (not C. coli or C. jejuni)
        - Enterococcus faecalis
        - Enterococcus faecium
        - Escherichia coli & Shigella spp.
        - Helicobacter pylori
        - Klebsiella spp.
        - Mycobacterium tuberculosis
        - Neisseria gonorrhoeae
        - Salmonella spp.
        - Staphylococcus aureus

    **XDR Shigella prediction**

    The `ResFinder` Task also has the ability to predict whether or not a sample meets the CDC's definition for extensively drug-resistant (XDR) Shigella. 

    > _CDC defines XDR Shigella bacteria as strains that are resistant to all commonly recommended empiric and alternative antibiotics — azithromycin, ciprofloxacin, ceftriaxone, trimethoprim-sulfamethoxazole (TMP-SMX), and ampicillin_. [See also the _Increase in Extensively Drug-Resistance Shigellosis in the United State_ CDC Health Network Alert](https://emergency.cdc.gov/han/2023/han00486.asp) where this definition can be found.
    
    ??? toggle "Criteria for XDR Shigella Prediction"
        A sample is required to meet **all 7 criteria** in order to be designated as `Potentially XDR Shigella` 

        1. The GAMBIT task in the workflow must identify the sample as `Shigella` OR the user must input the word `Shigella` somewhere within the input String variable `expected_taxon`. This requirement serves as the identification of a sample to be of the Shigella genus.
        2. Resfinder or PointFinder predicted resistance to **Ampicillin**
        3. Resfinder or PointFinder predicted resistance to **Azithromycin**
        4. Resfinder or PointFinder predicted resistance to **Ciprofloxacin**
        5. Resfinder or PointFinder predicted resistance to **Ceftriazone**
        6. Resfinder or PointFinder predicted resistance to **Trimethoprim**
        7. Resfinder or PointFinder predicted resistance to **Sulfamethoxazole**

    ??? toggle "Additional Criteria for XDR Shigella from CDC NARMS"

        CDC NARMS has a slightly different definition for XDR Shigella than the one used to determine the `resfinder_predicted_xdr_shigella` output above. Their definition is more stringent, and requires **two** different mechanisms that confer resistant to quinolones. The following two outputs can help determine if the sample meets this more strict criteria:

        - `resfinder_predicted_resistance_fq` - a semicolon-separated list of drugs (ciprofloxacin, fluoroquinolone, nalidixic acid, and unknown quinolone) and their associated mutations for which ResFinder and/or PointFinder predicted resistance.
        - `resfinder_predicted_resistance_fq_mechanisms` - the number of unique point mutations that confer fluoroquinolone resistance; mutations that confer resistance to more than one drug are counted once.

    There are 3 potential outputs for the **`resfinder_predicted_xdr_shigella`** output string:

    - **`Not Shigella based on gambit_predicted_taxon or user input`**
    - **`Not XDR Shigella`** for samples indicated as Shigella by either GAMBIT or user input _**BUT**_ ResFinder did **NOT** predict resistance to **all 6 drugs in the XDR definition above**
    - **`Potentially XDR Shigella`** for samples indicated as Shigella _**AND**_ ResFinder/PointFinder **DID** predict resistance to ceftriazone, azithromycin, ciprofloxacin, trimethoprim, sulfamethoxazole, _and_ ampicillin.
  
    !!! tip "What should I do if I see "Potentially XDR Shigella"?"
        If a sample is designated as `Potentially XDR Shigella`, we recommend reviewing the following output columns to confirm the presence of resistance genes/mutations associated with each drug in the XDR definition:
            
        - The output file found in `resfinder_pheno_table`
        - The output file found in `pointfinder_results` (if PointFinder was run)
        - The output string found in `resfinder_predicted_pheno_resistance`
    
        Feel free to contact us at <support@theiagen.com> for assistance in interpreting these results.
    

    !!! techdetails "ResFinder Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_resfinder.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/drug_resistance/task_resfinder.wdl) |
        | Software Source Code | [ResFinder Tool on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder/src/master/)<br>[ResFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)<br>[PointFinder Database on BitBucket](https://bitbucket.org/genomicepidemiology/pointfinder_db/src/master/) |
        | Software Documentation | [ResFinder on BitBucket](https://bitbucket.org/genomicepidemiology/resfinder/src/master/README.md) |
        | Original Publication(s) | _ResFinder tool_: [ResFinder 4.0 for predictions of phenotypes from genotypes](https://academic.oup.com/jac/article/75/12/3491/5890997)<br>_ResFinder database_: [Identification of acquired antimicrobial resistance genes](https://doi.org/10.1093/jac/dks261) |
