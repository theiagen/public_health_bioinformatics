??? task "`SeqSero2`: Serotyping"

    [SeqSero2](https://github.com/denglab/SeqSero2) is a tool for _Salmonella_ serotype prediction. In the TheiaProk Illumina workflows, SeqSero2 takes in raw sequencing reads and performs targeted assembly of serotype determinant alleles, which can be used to predict serotypes including contamination between serotypes. For the TheiaProk ONT and FASTA workflows, SeqSero2 uses the genome assembly as input.

    If reads are provided, SeqSero2 performs allele micro-assembly by default. This occurs through targeted assembly of serotype determinant alleles, and any assembled alleles are used to predict the sample's serotype, and can predict potential contamination. If the `seqsero2_mode` optional variable is changed to `"k"` (for k-mer mode), SeqSero2 will perform serotyping based on unique k-mers of serotype determinants. If the input data is an assembly FASTA, the k-mer mode must be used, and the genome assembly is used to generate the search k-mers instead of the raw reads. 
       
    !!! techdetails "SeqSero2 Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_seqsero2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_seqsero2.wdl) |
        | Software Source Code | [SeqSero2](https://github.com/denglab/SeqSero2) |
        | Software Documentation | [SeqSero2](https://github.com/denglab/SeqSero2) |
        | Original Publication(s) | [Salmonella serotype determination utilizing high-throughput genome sequencing data.](https://journals.asm.org/doi/10.1128/JCM.00323-15)<br>[SeqSero2: rapid and improved Salmonella serotype determination using whole genome sequencing data.](https://journals.asm.org/doi/10.1128/AEM.01746-19) |
