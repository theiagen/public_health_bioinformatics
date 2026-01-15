??? task "`quasitools`"

    `quasitools` performs genomic characterization for HIV by using the HyDRA module for identifying drug resistance mutations in HIV-1 samples based on the [Stanford HIV Drug Resistance Database](https://hivdb.stanford.edu/) and the [2009 WHO list for Surveillance of Transmitted HIVDR](https://hivdb.stanford.edu/page/who-sdrm-list/); see also the papers linked below.

    The HyDRA module in quasitools maps the sample sequence against an annotated HIV-1 reference and performs variant calling. Those variants are compared to the databases described above, and any matches are reported, along with the complete list of variants. 
    
    !!! techdetails "quasitools Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quasitools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/lentivirus/task_quasitools.wdl) |
        | Software Source Code | [quasitools on GitHub](https://github.com/phac-nml/quasitools/) |
        | Software Documentation | [quasitools HyDRA README](https://phac-nml.github.io/quasitools/hydra/) |
        | Original Publication(s) | quasitools preprint: _[quasitools: A Collection of Tools for Viral Quasispecies Analysis](https://doi.org/10.1101/733238)_<br>WHO 2009 Database: _[Drug resistance mutations for surveillance of transmitted HIV-1 drug-resistance: 2009 update](https://doi.org/10.1371/journal.pone.0004724)_<br>Stanford Database: _[Human immunodeficiency virus reverse transcriptase and protease sequence database](https://doi.org/10.1093/nar/gkg100)_ |
