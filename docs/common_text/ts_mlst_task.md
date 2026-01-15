??? task "`TS_MLST`: MLST Profiling"
    [Multilocus sequence typing (MLST)](https://doi.org/10.1073/pnas.95.6.3140) is a typing method reflecting population structure. It was developed as a portable, unambiguous method for global epidemiology using PCR, but can be applied to whole-genome sequences *in silico*. MLST is commonly used for pathogen surveillance, ruling out transmission, and grouping related genomes for comparative analysis.

    TheiaProk uses [the MLST tool developed by Torsten Seeman](https://github.com/tseemann/mlst) to assess MLST using traditional [PubMLST](https://pubmlst.org/) typing schemes. 

    MLST schemes are taxa-specific. Each scheme uses fragments of typically 7 housekeeping genes ("loci") and has a database associating an arbitrary number with each distinct allele of each locus. Each unique combination of alleles ("allelic profile") is assigned a numbered sequence type (ST). Significant diversification of genomes is captured by changes to the MLST loci via mutational events creating new alleles and STs, or recombinational events replacing the allele and changing the ST. Relationships between STs are based on the number of alleles they share. Clonal complexes share a scheme-specific number of alleles (usually for five of the seven loci).

    !!! tip "MLST Limitations"
        Some taxa have multiple MLST schemes, and some MLST schemes are insufficiently robust.
        
    ???+ toggle "Interpretation of MLST results"
        
        Each MLST results file returns the ST and allele results for one sample. If the alleles and ST are correctly assigned, only a single integer value will be present for each. If an ST cannot be assigned, multiple integers or additional characters will be shown, [representing the issues with assignment as described here](https://github.com/tseemann/mlst/tree/v2.22.0#missing-data).
        
    ??? toggle "Identifying novel alleles and STs"
        
        The MLST schemes used in TheiaProk are curated on the PubMLST website.If you identify novel alleles or allelic profiles in your data using TheiaProk's MLST task, you can get these assigned via PubMLST:
        
        1. Check that the novel allele or ST has not already been assigned a type on PubMLST. 
            1. Download the assembly file from Terra for your sample with the novel allele or ST
            2. Go to the [PubMLST webpage for the organism of interest](https://pubmlst.org/organisms) 
            3. Navigate to the organism "Typing" page
            4. Under "Query a sequence" choose "Single sequence" (e.g., [this is the page for _H. influenzae_](https://pubmlst.org/bigsdb?db=pubmlst_hinfluenzae_seqdef&page=sequenceQuery)), select the MLST scheme under "Please select locus/scheme", upload the assembly fasta file, and click submit
            5. Results will be returned lower on the page
        2. If the allele or ST has not been typed previously on the PubMLST website (step 1), new allele or ST numbers can be assigned using the instructions [provided by pubMLST](https://pubmlst.org/submit-data).
        
    ???+ toggle "Taxa with multiple MLST schemes"
        
        By default, the MLST tool automatically detects the genome's taxa to select the MLST scheme. Users may specify a desired scheme using the `scheme` variable of the `ts_mlst` task. Available schemes are listed [here](https://github.com/tseemann/mlst/blob/master/db/scheme_species_map.tab) and the scheme name should be provided in quotation marks ("...").
                
        Some taxa have multiple MLST schemes, e.g. the _Escherichia_ and _Leptospira_ genera,  _Acinetobacter baumannii_, _Clostridium difficile_, and _Streptococcus thermophilus_. Only one scheme will be used by default. This can be changed for _Escherichia_ and _A. baumannii_ by setting the `ts_mlst run_secondary_scheme` variable to `true`. This will cause MLST to run the other scheme associated with those organisms and output those results in addition to the default scheme.

        For _E. coli_, the user may set `ts_mlst scheme_override` to be `true` to prevent running _E. coli_ samples through the "aeromonas", "cfreundii", and "senterica" schemes, which can be common mischaracterizations.
                
    !!! techdetails "TS_MLST Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_ts_mlst.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/multi/task_ts_mlst.wdl) |
        | Software Source Code | [mlst on GitHub](https://github.com/tseemann/mlst) |
        | Software Documentation | [mlst on GitHub](https://github.com/tseemann/mlst) |
        | Original Publication(s) | _PubMLST database_: [Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications](https://doi.org/10.12688/wellcomeopenres.14826.1)
