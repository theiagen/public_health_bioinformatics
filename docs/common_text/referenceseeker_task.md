??? task "ReferenceSeeker Details (Optional)"
    ##### ReferenceSeeker
    `ReferenceSeeker` uses your draft assembly to identify closely related bacterial, viral, fungal, or plasmid genome assemblies in [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/).

    Databases that can be used with ReferenceSeeker are as follows, and can be used by pasting the GSURI in double quotation marks `" "` into the `referenceseeker_db` optional input:

      - archea:  `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-archaea-refseq-205.v20210406.tar.gz`
      - bacterial (**default**): `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-bacteria-refseq-205.v20210406.tar.gz`
      - fungi: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-fungi-refseq-205.v20210406.tar.gz`
      - plasmids: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-plasmids-refseq-205.v20210406.tar.gz`
      - viral: `gs://theiagen-public-files-rp/terra/theiaprok-files/referenceseeker-viral-refseq-205.v20210406.tar.gz`

    For ReferenceSeeker to identify a genome, it must meet user-specified thresholds for sequence coverage (`referenceseeker_conserved_dna_threshold`; default >= 0.69) and identity (`referenceseeker_ani_threshold`; default >= 0.95 ). 
    
    A list of closely related genomes is provided in `referenceseeker_tsv`. The reference genome that ranks highest according to ANI and conserved DNA values is considered the closest match and will be downloaded, with information about this provided in the `assembly_fetch_referenceseeker_top_hit_ncbi_accession` output.

    !!! techdetails "ReferenceSeeker Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_referenceseeker.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/phylogenetic_inference/utilities/task_referenceseeker.wdl) |
        | Software Source Code | [ReferenceSeeker on GitHub](https://github.com/oschwengers/referenceseeker) |
        | Software Documentation | [ReferenceSeeker on GitHub](https://github.com/oschwengers/referenceseeker) |
        | Original Publication(s) | [ReferenceSeeker: rapid determination of appropriate reference genomes](https://joss.theoj.org/papers/10.21105/joss.01994) |
