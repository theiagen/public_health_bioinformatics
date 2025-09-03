??? task "`LisSero`: Serogroup Prediction"
    [LisSero](https://github.com/MDU-PHL/LisSero) performs serogroup prediction for _Listeria monocytogenes_ based on the presence or absence of five genes, _lmo1118_, _lmo0737_, ORF2110, ORF2819, and _Prs_. LisSero does not predict somatic (O) or flagellar (H) biosynthesis.

    LisSero uses BLAST to query the sample against a built-in database to determine the most likely serogroup. 

    ??? toggle "Serogroup Decision Tree"
        The following decision tree is used to determine the serogroup for the sample.
        ```
        if not Prs:
            stop
        elif lmo0737 and not (ORF2819 or ORF2110):
            if lmo1118:
                Serogroup 1/2c, 3c
            else:
                Serogroup 1/2a, 3a
        elif ORF2819 and not (lmo0737 or lmo1118):
            if ORF2110:
                Serogroup 4b, 4d, 4e
            else:
                Serogroup 1/2b, 3b, 7
        elif lmo0737 and ORF2819 and ORF2110:
            Serogroup 4b, 4d, 4e*
        else:
            Nontypable
        ```

    !!! techdetails "LisSero Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_lissero.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/listeria/task_lissero.wdl) |
        | Software Source Code | [LisSero](https://github.com/MDU-PHL/LisSero) |
        | Software Documentation | [LisSero](https://github.com/MDU-PHL/LisSero) |
