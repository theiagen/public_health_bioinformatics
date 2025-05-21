??? task "`Flye`: _De novo_ Assembly"
    `flye_denovo` is a sub-workflow that performs _de novo_ assembly using Flye for ONT data and supports additional polishing and visualization steps.
    
    !!! tip "Ensure correct medaka model is selected if performing medaka polishing"
        In order to obtain the best results, the appropriate model must be set to match the sequencer's basecaller model; this string takes the format of {pore}\_{device}\_{caller variant}\_{caller_version}. See also <https://github.com/nanoporetech/medaka?tab=readme-ov-file#models>. If `flye` is being run on legacy data the medaka model will likely be `r941_min_hac_g507`. Recently generated data will likely be suited by the default model of `r1041_e82_400bps_sup_v5.0.0`.

    The detailed steps and tasks are as follows:

    ??? toggle "`Porechop`: Read Trimming (optional; off by default)"
        Read trimming is optional and can be enabled by setting the `run_porchop` input variable to true.

        Porechop is a tool for finding and removing adapters from ONT data. Adapters on the ends of reads are trimmed, and when a read has an adapter in the middle, the read is split into two.

        !!! techdetails "Porechop Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_porechop.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_porechop.wdl) |
            | Software Source Code | [Porechop on GitHub](https://github.com/rrwick/Porechop) |
            | Software Documentation | [https://github.com/rrwick/Porechop#porechop](https://github.com/rrwick/Porechop#porechop) |

{{ include_md("common_text/flye_task.md", indent=4, replacements={'??? task "`flye`"' : '??? toggle "`Flye`: _De novo_ Assembly"'}) }}

    ??? toggle "`Bandage`: Graph Visualization"
        Bandage creates _de novo_ assembly graphs containing the assembled contigs and the connections between those contigs.

        !!! techdetails "Bandage Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_bandage_plot.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_bandage_plot.wdl) |
            | Software Source Code | [Bandage on GitHub](https://github.com/rrwick/Bandage) |
            | Software Documentation | [Bandage Documentation](https://github.com/rrwick/Bandage#bandage) |
            | Original Publication(s) | [Bandage: interactive visualization of _de novo_ genome assemblies](https://academic.oup.com/bioinformatics/article/31/20/3350/196114) |

    ??? toggle "`Polypolish`: Hybrid Assembly Polishing ==_for ONT and Illumina data_=="
        If short reads are provided with the optional `illumina_read1` and `illumina_read2` inputs, Polypolish will use those short-reads to correct errors in the long-read assemblies. Uniquely, Polypolish uses the short-read alignments where each read is aligned to _all_ possible locations, meaning that repeat regions will have error correction.
    
        !!! techdetails "Polypolish Technical Details"
            |  | Links |
            | --- | --- |
            | Task | [task_polypolish.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_polypolish.wdl) |
            | Software Source Code | [Polypolish on GitHub](https://github.com/rrwick/Polypolish) |
            | Software Documentation | [Polypolish Documentation](https://github.com/rrwick/Polypolish#polypolish) |
            | Original Publication(s) | [Polypolish: short-read polishing of long-read bacterial genome assemblies](https://doi.org/10.1371/journal.pcbi.1009802)<br>[How low can you go? Short-read polishing of Oxford Nanopore bacterial genome assemblies](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.001254) |

    ??? toggle "`Medaka`: Polishing of Flye assembly (default; optional)"
        Polishing is optional and can be skipped by setting the `skip_polishing` variable to true. If polishing is skipped, then neither Medaka or Racon will run.

        Medaka is the default assembly polisher used in TheiaProk. Racon may be used alternatively, and if so, Medaka will not run. Medaka uses the raw reads to polish the assembly and generate a consensus sequence. 

        Importantly, Medaka requires knowing the model that was used to generate the read data. There are several ways to provide this information:

        - Automatic Model Selection: Automatically determines the most appropriate Medaka model based on the input data, ensuring optimal polishing results without manual intervention. 
        - User-Specified Model Override: Allows users to specify a particular `Medaka model` if automatic selection does not yield the desired outcome or for specialized use cases.
        - Default Model: If both automatic model selection fails and no user-specified model is provided, Medaka defaults to the predefined fallback model `r1041_e82_400bps_sup_v5.0.0`. 

        !!! info "Medaka Model Resolution Process" 
            Medaka's automatic model selection uses the `medaka tools resolve_model` command to identify the appropriate model for polishing. This process relies on metadata embedded in the input file, which is typically generated by the basecaller. If the automatic selection fails to identify a suitable model, Medaka gracefully falls back to the default model to maintain workflow continuity. **Users should verify the chosen model and consider specifying a model override if necessary.**

        !!! techdetails "Medaka Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_medaka.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_medaka.wdl) |
            | Software Source Code | [Medaka on GitHub](https://github.com/nanoporetech/medaka) |
            | Software Documentation | [Medaka Documentation](https://github.com/nanoporetech/medaka#medaka) |

    ??? toggle "`Racon`: Polishing of Flye assembly (alternative; optional)"
        Polishing is optional and can be skipped by setting the `skip_polishing` variable to true. If polishing is skipped, then neither Medaka or Racon will run.

        `Racon` is an alternative to using `medaka` for assembly polishing, and can be run by setting the `polisher` input to "racon".  Racon is a consensus algorithm designed for refining raw de novo DNA assemblies generated from long, uncorrected sequencing reads.

        !!! techdetails "Racon Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_racon.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/polishing/task_racon.wdl) |
            | Software Source Code | [Racon on GitHub](https://github.com/lbcb-sci/racon) |
            | Software Documentation | [Racon Documentation](https://github.com/lbcb-sci/racon#racon) |
            | Original Publication(s) | [Fast and accurate de novo genome assembly from long uncorrected reads](https://genome.cshlp.org/content/27/5/737) |

    ??? toggle "`Filter Contigs`: Filter contigs below a threshold length and remove homopolymer contigs"
        This task filters the created contigs based on a user-defined minimum length threshold (default of 1000) and eliminates homopolymer contigs (contigs of any length that consist of a single nucleotide). This ensures high-quality assemblies by retaining only contigs that meet specified criteria. Detailed metrics on contig counts and sequence lengths before and after filtering are provided in the output.

        !!! techdetails "Filter Contigs Technical Details" 
            |  | Links |
            | --- | --- |
            | WDL Task | [task_filter_contigs.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_filter_contigs.wdl) |
        
    ??? toggle "`Dnaapler`: Final Assembly Orientation"
        Dnaapler reorients contigs to start at specific reference points. Dnaapler supports the following modes, which can be indicated by filling the `dnaapler_mode` input variable with the desired mode. The default is `all`, which reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474`.

        - **all**: Reorients contigs to start with `dnaA`, `terL`, `repA`, or `COG1474` (_Default_)
        - **chromosome**: Reorients to begin with the `dnaA` chromosomal replication initiator gene, commonly used for bacterial chromosome assemblies.
        - **plasmid**: Reorients to start with the `repA` plasmid replication initiation gene, ideal for plasmid assemblie
        - **phage**: Reorients to start with the `terL` large terminase subunit gene, used for bacteriophage assemblies
        - **archaea**: Reorients to start with the `COG1474` archaeal Orc1/cdc6 gene, relevant for archaeal assemblies
        - **custom**: Reorients based on a user-specified gene in amino acid FASTA format for experimental or unique workflows
        - **mystery**: Reorients to start with a random CDS for exploratory purposes
        - **largest**: Reorients to start with the largest CDS in the assembly, often useful for poorly annotated genomes
        - **nearest**: Reorients to start with the first CDS nearest to the sequence start, resolving CDS breakpoints
        - **bulk**: Processes multiple contigs to start with the desired start gene (`dnaA`, `terL`, `repA`, or custom)

        !!! techdetails "Dnaapler Technical Details"
            |  | Links |
            | --- | --- |
            | WDL Task | [task_dnaapler.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_dnaapler.wdl) |
            | Software Source Code | [Dnaapler on GitHub](https://github.com/gbouras13/dnaapler) |
            | Software Documentation | [Dnaapler Documentation](https://github.com/gbouras13/dnaapler?tab=readme-ov-file#dnaapler) |
            | Original Publication(s) | [Dnaapler: a tool to reorient circular microbial genomes](https://joss.theoj.org/papers/10.21105/joss.05968) |