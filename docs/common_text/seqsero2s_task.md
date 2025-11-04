??? task "`SeqSero2S`: Serotyping"
    [SeqSero2S](https://github.com/LSTUGA/SeqSero2S) is a tool for _Salmonella_ serotype prediction. SeqSero2S is a software package that determines serotype antigens by finding the genes responsible for the flagellar filament (H antigen) and the outermost oligosaccharides of LPS (O antigen) in _Salmonella_ and matches them to known representative alleles for those genes in a curated database. An antigenic formula describes all serotypes; _Salmonella enterica_ subspecies _enterica_ serotypes are also assigned a serotype name by the genotyping tool.

    In the TheiaProkIllumina workflows, SeqSero2S takes in raw sequencing reads and performs targeted assembly of serotype determinant alleles, which can be used to predict serotypes including contamination between serotypes. For the TheiaProk ONT and FASTA workflows, SeqSero2S uses the genome assembly as input.

    If reads are provided, SeqSero2S performs allele micro-assembly by default. This occurs through targeted assembly of serotype determinant alleles, and any assembled alleles are used to predict the sample's serotype, and can predict potential contamination. If the `seqsero2s_mode` optional variable is changed to `"k"` (for k-mer mode), SeqSero2S will perform serotyping based on unique k-mers of serotype determinants. If the input data is an assembly FASTA, the k-mer mode _must_ be used, and the genome assembly is used to generate the search k-mers instead of the raw reads. 

    ??? dna "What is the difference between SeqSero2 and SeqSero2S?"
        Recently, the SeqSero2S software was released and adopted by the PulseNet 2.0 (PN2.0) platform for subspecies identification and serotype determination of _Salmonella_ spp. Genetic determination of rarer serotypes can be problematic due to a lack of sequences for rare antigen types and alleles, a lack of understanding of the genetic basis for some antigens, or some inconsistencies in the White-Kauffmann-Le Minor (WKL) Scheme for _Salmonella_ serotype designation. As such, SeqSero2S predicts serotypes using a simplified interpretation based on the most commonly seen serotypes. SeqSero2S also includes additional functionalities to improve genomic prediction of serotypes in general and mitigate potential drawbacks of the simplified scheme.
        
        The simplification of SeqSero2S can be summarized as such:

        - 178 serotypes were provisionally removed
            - 84 serotypes, of which the exact antigenic formula could not be determined
            - 66 serotypes had no alleles, or probes were not available
            - 10 serotypes had no probes available for O:9,46,27
            - 14 serotypes were removed because O:54 is encoded on a plasmid
            - 4 serotypes could not be differentiated from other serotypes
        - 57 serotypes were merged into other serotypes
        - The antigenic formulae of 890 serotypes were simplified, but the simplifications did not result in merging
       
    !!! techdetails "SeqSero2S Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_seqsero2s.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/salmonella/task_seqsero2s.wdl) |
        | Software Source Code | [SeqSero2S on GitHub](https://github.com/LSTUGA/SeqSero2S) |
        | Software Documentation | [SeqSero2S on GitHub](https://github.com/LSTUGA/SeqSero2S) |
        | Original Publication(s) | [Salmonella serotype determination utilizing high-throughput genome sequencing data.](https://journals.asm.org/doi/10.1128/JCM.00323-15)<br>[SeqSero2: rapid and improved Salmonella serotype determination using whole genome sequencing data.](https://journals.asm.org/doi/10.1128/AEM.01746-19) |
