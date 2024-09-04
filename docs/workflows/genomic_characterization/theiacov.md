# TheiaCoV Workflow Series

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line compatibliity** | **Workflow type** |
|---|---|---|---|---|
| [Genomic Characterization](../../workflows_overview/workflows-type.md/#genomic-characterization) | [Viral](../../workflows_overview/workflows-kingdom.md/#viral) | PHB v2.2.0 | Yes | Sample-level |

## TheiaCoV Workflows

!!! dna inline end "Key Resources"

    [**Reference Materials for SARS-CoV-2**](https://www.notion.so/Docker-Image-and-Reference-Materials-for-SARS-CoV-2-Genomic-Characterization-98328c61f5cb4f77975f512b55d09108?pvs=21)

    [**Reference Materials for Mpox**](https://www.notion.so/Workspace-Reference-Materials-for-MPXV-Genomic-Characterization-a34f355c68c54c0a82e926d4de607bca?pvs=21)

    ??? toggle "HIV Input JSONs"
        - [TheiaCoV_Illumina_PE_HIV_v1_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/439f1c74-d91e-4978-b173-3302f878e343/TheiaCoV_Illumina_PE_HIV_v1_2024-04-19.json)
        - [TheiaCoV_Illumina_PE_HIV_v2_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/2c7872de-44c8-406d-bbec-fadaacbb0d98/TheiaCoV_Illumina_PE_HIV_v2_2024-04-19.json)
        - [TheiaCoV_ONT_HIV_v1_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/9f9a7bd1-2ac4-47fb-967b-4198a45d4a71/TheiaCoV_ONT_HIV_v1_2024-04-19.json)
        - [TheiaCoV_ONT_HIV_v2_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/13fdfec0-4a81-460e-948a-be6ad30d022d/TheiaCoV_ONT_HIV_v2_2024-04-19.json)
        
    ??? toggle "WNV Input JSONs"
        - [TheiaCoV_Illumina_PE_WNV_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/6af74d02-9985-428d-897e-e04ebacc42a3/TheiaCoV_Illumina_PE_WNV_2024-04-19.json)
        - [TheiaCoV_Illumina_SE_WNV_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/cb8dec19-2563-4070-9ae9-031c089f8b3d/TheiaCoV_Illumina_SE_WNV_2024-04-19.json)
        - [TheiaCoV_FASTA_WNV_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/f2059069-5ce1-45e1-ab9e-51925158c0eb/TheiaCoV_FASTA_WNV_2024-04-19.json)
        
    ??? toggle "Flu Input JSONs"
        - [TheiaCoV_Illumina_PE_flu_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/ba326b69-8a2a-4af2-a74f-e710e667f82b/TheiaCoV_Illumina_PE_flu_2024-04-19.json)
        - [TheiaCoV_ONT_flu_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/c01c98f5-d00e-4ff2-ad09-6cc3ff1ad3a7/TheiaCoV_ONT_flu_2024-04-19.json)
        - [TheiaCoV_FASTA_flu_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/4c7d7a16-2c20-4cbc-9618-231afade9940/TheiaCoV_FASTA_flu_2024-04-19.json)
        
    ??? toggle "RSV-A Input JSONs"
        - [TheiaCoV_Illumina_PE_RSV-B_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/2be20bb8-b733-4f02-a27f-b0cf19d015f8/TheiaCoV_Illumina_PE_RSV-B_2024-04-19.json)
        - [TheiaCoV_FASTA_RSV-A_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/ba6a4845-14ee-4664-b9f3-808f76c87d15/TheiaCoV_FASTA_RSV-A_2024-04-19.json)
        
    ??? toggle "RSV-B Input JSONs" 
        - [TheiaCoV_Illumina_PE_RSV-A_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/dd1612ff-20c5-4310-9cb3-c07bf9b7e8a1/TheiaCoV_Illumina_PE_RSV-A_2024-04-19.json)
        - [TheiaCoV_FASTA_RSV-B_2024-04-19.json](https://prod-files-secure.s3.us-west-2.amazonaws.com/be290196-9090-4f3c-a9ab-fe730ad213e0/160cdfbc-a556-40bc-aa05-84ae69511400/TheiaCoV_FASTA_RSV-B_2024-04-19.json)

**The TheiaCoV workflows are for the assembly, quality assessment, and characterization of viral genomes.** There are currently five TheiaCoV workflows designed to accommodate different kinds of input data:

1. Illumina paired-end sequencing (**TheiaCoV_Illumina_PE**)
2. Illumina single-end sequencing (**TheiaCoV_Illumina_SE**)
3. ONT sequencing (**TheiaCoV_ONT**)
4. Genome assemblies (**TheiaCoV_FASTA**)
5. ClearLabs sequencing (**TheiaCoV_ClearLabs**)

Additionally, the **TheiaCoV_FASTA_Batch** workflow is available to process several hundred SARS-CoV-2 assemblies at the same time.

---

### Supported Organisms

These workflows currently support the following organisms:

- **SARS-CoV-2** (`"sars-cov-2"`, `"SARS-CoV-2"`) - ==_default organism input_==
- **Monkeypox virus** (`"MPXV"`, `"mpox"`, `"monkeypox"`, `"Monkeypox virus"`, `"Mpox"`)
- **Human Immunodeficiency Virus** (`"HIV"`)
- **West Nile Virus** (`"WNV"`, `"wnv"`, `"West Nile virus"`)
- **Influenza** (`"flu"`, `"influenza"`, `"Flu"`, `"Influenza"`)
- **RSV-A** (`"rsv_a"`, `"rsv-a"`, `"RSV-A"`, `"RSV_A"`)
- **RSV-B** (`"rsv_b"`, `"rsv-b"`, `"RSV-B"`, `"RSV_B"`)

The compatibility of each workflow with each pathogen is shown below:

|  | SARS-CoV-2 | Mpox | HIV | WNV | Influenza | RSV-A | RSV-B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Illumina_PE | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| Illumina_SE | ✅ | ✅ | ❌ | ✅ | ❌ | ✅ | ✅ |
| ClearLabs | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| ONT | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| FASTA | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ |

We’ve provided the following information to help you set up the workflow for each organism in the form of input JSONs.

### Inputs

[TheiaCoV Inputs](https://www.notion.so/cfe814ea451e404bb11af6be2333b97a?pvs=21)

[TheiaCoV_FASTA_Batch Inputs](https://www.notion.so/3316ffb500524eafaa0b0b9167bab919?pvs=21)

!!! dna ""
    ??? toggle "TheiaCoV_Illumina_PE Input Read Data"

        The TheiaCoV_Illumina_PE workflow takes in Illumina paired-end read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before Terra uploads to minimize data upload time.

        By default, the workflow anticipates **2 x 150bp** reads (i.e. the input reads were generated using a 300-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate shorter read data, such as the 2 x 75bp reads generated using a 150-cycle sequencing kit.

    ??? toggle "TheiaCoV_Illumina_SE Input Read Data"

        TheiaCoV_Illumina_SE takes in Illumina single-end reads. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. Theiagen highly recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time & save on storage costs.

        By default, the workflow anticipates **1 x 35 bp** reads  (i.e. the input reads were generated using a 70-cycle sequencing kit). Modifications to the optional parameter for `trim_minlen` may be required to accommodate longer read data.

    ??? toggle "TheiaCoV_ONT Input Read Data"

        The TheiaCoV_ONT workflow takes in base-called ONT read data. Read file names should end with `.fastq` or `.fq`, with the optional addition of `.gz`. When possible, Theiagen recommends zipping files with [gzip](https://www.gnu.org/software/gzip/) before uploading to Terra to minimize data upload time.

        **The ONT sequencing kit and base-calling approach can produce substantial variability in the amount and quality of read data. Genome assemblies produced by the TheiaCoV_ONT workflow must be quality assessed before reporting results.**

    ??? toggle "TheiaCoV_FASTA Input Assembly Data"

        The TheiaCoV_FASTA workflow takes in assembly files in FASTA format.

    ??? toggle "TheiaCoV_ClearLabs Input Read Data"

        The TheiaCoV_ClearLabs workflow takes in read data produced by the Clear Dx platform from ClearLabs. However, many users use the TheiaCoV_FASTA workflow instead of this one due to a few known issues when generating assemblies with this pipeline that are not present when using ClearLabs-generated FASTA files.

    ??? toggle "TheiaCoV_FASTA_Batch Input Data"

        The TheiaCoV_FASTA_Batch workflow takes in a set of assembly files in FASTA format.

### Organism-specific parameters and logic {#org-specific}

The `organism_parameters` sub-workflow is the first step in all TheiaCoV workflows. This step automatically sets the different parameters needed for each downstream tool to the appropriate value for the user-designated organism (by default, `"sars-cov-2"` is the default organism).

!!! dna ""
    The following tables include the relevant organism-specific parameters; **all of these default values can be overwritten by providing a value for the "Overwrite Variable Name" field**.

    ??? toggle "SARS-CoV-2 Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** |
        |---|---|---|
        | gene_locations_bed_file | sars-cov-2 | `"gs://theiagen-public-files-rp/terra/sars-cov-2-files/sc2_gene_locations.bed"` |
        | genome_length_input | sars-cov-2 | `29903` |
        | nextclade_dataset_name_input | sars-cov-2 | `"nextstrain/sars-cov-2/wuhan-hu-1/orfs"` |
        | nextclade_dataset_tag_input | sars-cov-2 | `"2024-07-17--12-57-03Z"` |
        | pangolin_docker_image | sars-cov-2 | `"us-docker.pkg.dev/general-theiagen/staphb/pangolin:4.3.1-pdata-1.29 "`|
        | reference_genome | sars-cov-2 | `"gs://theiagen-public-files-rp/terra/augur-sars-cov-2-references/MN908947.fasta"` |
        | vadr_max_length | sars-cov-2 | `30000` |
        | vadr_mem | sars-cov-2 | `8` |
        | vadr_options | sars-cov-2 | `"--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"` |

    ??? toggle "Mpox Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** |
        |---|---|---|
        | gene_locations_bed_file | MPXV | `"gs://theiagen-public-files/terra/mpxv-files/mpox_gene_locations.bed"` |
        | genome_length_input | MPXV | `197200` |
        | kraken_target_organism_input | MPXV | `"Monkeypox virus"` |
        | nextclade_dataset_name_input | MPXV | `"nextstrain/mpox/lineage-b.1"` |
        | nextclade_dataset_tag_input | MPXV | `"2024-04-19--07-50-39Z"` |
        | primer_bed_file | MPXV | `"gs://theiagen-public-files/terra/mpxv-files/MPXV.primer.bed"` |
        | reference_genome | MPXV | `"gs://theiagen-public-files/terra/mpxv-files/MPXV.MT903345.reference.fasta"` |
        | reference_gff_file | MPXV | `"gs://theiagen-public-files/terra/mpxv-files/Mpox-MT903345.1.reference.gff3"` |
        | vadr_max_length | MPXV | `210000` |
        | vadr_mem | MPXV | `8` |
        | vadr_options | MPXV | `"--glsearch -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --out_allfasta --minimap2 --s_overhang 150"` |

    ??? toggle "WNV Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | genome_length_input | WNV | `11000` | |
        | kraken_target_organism_input | WNV | `"West Nile virus`" | |
        | nextclade_dataset_name_input | WNV | `"NA"` | TheiaCoV's Nextclade currently does not support WNV |
        | nextclade_dataset_tag_input | WNV | `"NA"` | TheiaCoV's Nextclade currently does not support WNV |
        | primer_bed_file | WNV | `"gs://theiagen-public-files/terra/theiacov-files/WNV/WNV-L1_primer.bed"` |  |
        | reference_genome | WNV | `"gs://theiagen-public-files/terra/theiacov-files/WNV/NC_009942.1_wnv_L1.fasta"` |  |
        | vadr_max_length | WNV | `11000` |  |
        | vadr_mem | WNV | `8` |  |
        | vadr_options | WNV | `"--mkey flavi --mdir /opt/vadr/vadr-models-flavi/ --nomisc --noprotid --out_allfasta"` |  |

    ??? toggle "Flu Defaults"
        | **Overwrite Variable Name** | **Organism** | **Flu Segment** | **Flu Subtype** | **Default Value** | **Notes** |
        |---|---|---|---|---|---|
        | flu_segment | flu | all | all | N/A | TheiaCoV will attempt to automatically assign a flu segment  |
        | flu_subtype | flu | all | all | N/A | TheiaCoV will attempt to automatically assign a flu subtype |
        | genome_length_input | flu | all | all | `13500` |  |
        | vadr_max_length | flu | all | all | `13500` |  |
        | vadr_mem | flu | all | all | `8` |  |
        | vadr_options | flu | all | all | `"--atgonly --xnocomp --nomisc --alt_fail extrant5,extrant3 --mkey flu"` |  |
        | nextclade_dataset_name_input | flu | ha | h1n1 | `"nextstrain/flu/h1n1pdm/ha/MW626062"` |  |
        | nextclade_dataset_tag_input | flu | ha | h1n1 | `"2024-07-03--08-29-55Z"` |  |
        | reference_genome | flu | ha | h1n1 | `"gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_ha.fasta"` |  |
        | nextclade_dataset_name_input | flu | ha | h3n2 | `"nextstrain/flu/h3n2/ha/EPI1857216"` |  |
        | nextclade_dataset_tag_input | flu | ha | h3n2 | `"2024-08-08--05-08-21Z"` |  |
        | reference_genome | flu | ha | h3n2 | `"gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_ha.fasta"` |  |
        | nextclade_dataset_name_input | flu | ha | victoria | `"nextstrain/flu/vic/ha/KX058884"` |  |
        | nextclade_dataset_tag_input | flu | ha | victoria | `"2024-07-03--08-29-55Z"` |  |
        | reference_genome | flu | ha | victoria | `"gs://theiagen-public-files-rp/terra/flu-references/reference_vic_ha.fasta"` |  |
        | nextclade_dataset_name_input | flu | ha | yamagata | `"nextstrain/flu/yam/ha/JN993010"` |  |
        | nextclade_dataset_tag_input | flu | ha | yamagata | `"2024-01-30--16-34-55Z"` |  |
        | reference_genome | flu | ha | yamagata | `"gs://theiagen-public-files-rp/terra/flu-references/reference_yam_ha.fasta"` |  |
        | nextclade_dataset_name_input | flu | na | h1n1 | `"nextstrain/flu/h1n1pdm/na/MW626056"` |  |
        | nextclade_dataset_tag_input | flu | na | h1n1 | `"2024-07-03--08-29-55Z"` |  |
        | reference_genome | flu | na | h1n1 | `"gs://theiagen-public-files-rp/terra/flu-references/reference_h1n1pdm_na.fasta"` |  |
        | nextclade_dataset_name_input | flu | na | h3n2 | `"nextstrain/flu/h3n2/na/EPI1857215"` |  |
        | nextclade_dataset_tag_input | flu | na | h3n2 | `"2024-04-19--07-50-39Z"` |  |
        | reference_genome | flu | na | h3n2 | `"gs://theiagen-public-files-rp/terra/flu-references/reference_h3n2_na.fasta"` |  |
        | nextclade_dataset_name_input | flu | na | victoria | `"nextstrain/flu/vic/na/CY073894"` |  |
        | nextclade_dataset_tag_input | flu | na | victoria | `"2024-04-19--07-50-39Z"` |  |
        | reference_genome | flu | na | victoria | `"gs://theiagen-public-files-rp/terra/flu-references/reference_vic_na.fasta"` |  |
        | nextclade_dataset_name_input | flu | na | yamagata | `"NA"` |  |
        | nextclade_dataset_tag_input | flu | na | yamagata | `"NA"` |  |
        | reference_genome | flu | na | yamagata | `"gs://theiagen-public-files-rp/terra/flu-references/reference_yam_na.fasta"` |  |

    ??? toggle "RSV-A Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | genome_length_input | rsv_a | 16000 |  |
        | kraken_target_organism | rsv_a | Respiratory syncytial virus |  |
        | nextclade_dataset_name_input | rsv_a | nextstrain/rsv/a/EPI_ISL_412866 |  |
        | nextclade_dataset_tag_input | rsv_a | 2024-08-01--22-31-31Z |  |
        | reference_genome | rsv_a | gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_a.fasta |  |
        | vadr_max_length | rsv_a | 15500 |  |
        | vadr_mem | rsv_a | 32 |  |
        | vadr_options | rsv_a | -r --mkey rsv --xnocomp |  |

    ??? toggle "RSV-B Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | genome_length_input | rsv_b | 16000 |  |
        | kraken_target_organism | rsv_b |  "Human orthopneumovirus" |  |
        | nextclade_dataset_name_input | rsv_b | nextstrain/rsv/b/EPI_ISL_1653999 | |
        | nextclade_dataset_tag_input | rsv_b | "2024-08-01--22-31-31Z" |  |
        | reference_genome | rsv_b | gs://theiagen-public-files-rp/terra/rsv_references/reference_rsv_b.fasta |  |
        | vadr_max_length | rsv_b | 15500 |  |
        | vadr_mem | rsv_b | 32 |  |
        | vadr_options | rsv_b | -r --mkey rsv --xnocomp |  |

    ??? toggle "HIV Defaults"
        | **Overwrite Variable Name** | **Organism** | **Default Value** | **Notes** |
        |---|---|---|---|
        | kraken_target_organism_input | HIV | Human immunodeficiency virus 1 |  |
        | genome_length_input | HIV-v1 | 9181 | This version of HIV originates from Oregon |
        | primer_bed_file | HIV-v1 | gs://theiagen-public-files/terra/hivgc-files/HIV-1_v1.0.primer.hyphen.bed | This version of HIV originates from Oregon |
        | reference_genome | HIV-v1 | gs://theiagen-public-files/terra/hivgc-files/NC_001802.1.fasta | This version of HIV originates from Oregon |
        | reference_gff_file | HIV-v1 | gs://theiagen-public-files/terra/hivgc-files/NC_001802.1.gff3 | This version of HIV originates from Oregon |
        | genome_length_input | HIV-v2 | 9840 | This version of HIV originates from Southern Africa |
        | primer_bed_file | HIV-v2 | gs://theiagen-public-files/terra/hivgc-files/HIV-1_v2.0.primer.hyphen400.1.bed | This version of HIV originates from Southern Africa |
        | reference_genome | HIV-v2 | gs://theiagen-public-files/terra/hivgc-files/AY228557.1.headerchanged.fasta | This version of HIV originates from Southern Africa |
        | reference_gff_file | HIV-v2 | gs://theiagen-public-files/terra/hivgc-files/AY228557.1.gff3 | This version of HIV originates from Southern Africa |

### Workflow Tasks

All input reads are processed through "core tasks" in the TheiaCoV Illumina, ONT, and ClearLabs workflows. These undertake read trimming and assembly appropriate to the input data type. TheiaCoV workflows subsequently launch default genome characterization modules for quality assessment, and additional taxa-specific characterization steps. When setting up the workflow, users may choose to use "optional tasks" as additions or alternatives to tasks run in the workflow by default.

#### Core tasks

!!! tip ""
    These tasks are performed regardless of organism, and perform read trimming and various quality control steps.

??? task "`versioning`: Version capture for TheiaEuk"

    The `versioning` task captures the workflow version from the GitHub (code repository) version.
        
    !!! techdetails "Version Capture Technical details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_versioning.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/task_versioning.wdl) |

??? task "`screen`: Total Raw Read Quantification and Genome Size Estimation"

    The [`screen`](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_screen.wdl) task ensures the quantity of sequence data is sufficient to undertake genomic analysis. It uses bash commands for quantification of reads and base pairs, and [mash](https://mash.readthedocs.io/en/latest/index.html) sketching to estimate the genome size and its coverage. At each step, the results are assessed relative to pass/fail criteria and thresholds that may be defined by optional user inputs. Samples that do not meet these criteria will not be processed further by the workflow:

    1. Total number of reads: A sample will fail the read screening task if its total number of reads is less than or equal to `min_reads`.
    2. The proportion of basepairs reads in the forward and reverse read files: A sample will fail the read screening if fewer than `min_proportion` basepairs are in either the reads1 or read2 files.
    3. Number of basepairs: A sample will fail the read screening if there are fewer than `min_basepairs` basepairs
    4. Estimated genome size:  A sample will fail the read screening if the estimated genome size is smaller than `min_genome_size` or bigger than `max_genome_size`.
    5. Estimated genome coverage: A sample will fail the read screening if the estimated genome coverage is less than the `min_coverage`.

    Read screening is undertaken on both the raw and cleaned reads. The task may be skipped by setting the `skip_screen` variable to true.

    Default values vary between the PE and SE workflow. The rationale for these default values can be found below.
    
    | Variable  | Rationale |
    | --- | --- |
    | `skip_screen` | Prevent the read screen from running |
        | `skip_screen` | Saving waste of compute resources on insufficient data |
        | `min_reads` | Minimum number of base pairs for 10x coverage of the Hepatitis delta (of the *Deltavirus* genus) virus divided by 300 (longest Illumina read length) |
        | `min_basepairs` | Greater than 10x coverage of the Hepatitis delta (of the *Deltavirus* genus) virus |
        | `min_genome_size` | Based on the Hepatitis delta (of the *Deltavirus* genus) genome- the smallest viral genome as of 2024-04-11 (1,700 bp) |
        | `max_genome_size` | Based on the *Pandoravirus salinus* genome, the biggest viral genome, (2,673,870 bp) with 2 Mbp added |
        | `min_coverage` | A bare-minimum coverage for genome characterization. Higher coverage would be required for high-quality phylogenetics. |
        | `min_proportion` | Greater than 50% reads are in the read1 file; others are in the read2 file |

    !!! techdetails "Screen Technical Details"
        
        There is a single WDL task for read screening. The `screen` task is run twice, once for raw reads and once for clean reads.
        
        |  | Links |
        | --- | --- |
        | Task | [task_screen.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_screen.wdl)  |

??? task "`read_QC_trim_pe` and `read_QC_trim_se`: Read Quality Trimming, Host and Adapter Removal, Quantification, and Identification _for Illumina workflows_"

    `read_QC_trim` is a sub-workflow within TheiaCoV that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below. The differences between TheiaCoV PE and SE in the `read_QC_trim` sub-workflow lie in the default parameters, the use of two or one input read file(s), and the different output files.

    ??? toggle "Host removal"
    
        All reads of human origin **are removed**, including their mates, by using NCBI’s [**human read removal tool (HRRT)**](https://github.com/ncbi/sra-human-scrubber). 

        HRRT is based on the [SRA Taxonomy Analysis Tool](https://doi.org/10.1186/s13059-021-02490-0) and employs a k-mer database constructed of k-mers from Eukaryota derived from all human RefSeq records with any k-mers found in non-Eukaryota RefSeq records subtracted from the database.

        !!! techdetails "NCBI-Scrub Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_ncbi_scrub.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_ncbi_scrub.wdl) |
            | Software source code | [NCBI Scrub on GitHub](https://github.com/ncbi/sra-human-scrubber) |
            | Software documentation | <https://github.com/ncbi/sra-human-scrubber/blob/master/README.md> |

    ??? toggle "Read quality trimming"

        Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_minlen`. 

        If fastp is selected for analysis, fastp also implements the additional read-trimming steps indicated below:

        | **Parameter** | **Explanation** |
        | --- | --- |
        | -g | enables polyG tail trimming |
        | -5 20 | enables read end-trimming |
        | -3 20 | enables read end-trimming |
        | --detect_adapter_for_pe | enables adapter-trimming **only for paired-end reads** |

    ??? toggle "Adapter removal"

        The `BBDuk` task removes adapters from sequence reads. To do this:

        - [Repair](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/) from the [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) package reorders reads in paired fastq files to ensure the forward and reverse reads of a pair are in the same position in the two fastq files.
        - [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)  (*"Bestus Bioinformaticus" Decontamination Using Kmers*) is then used to trim the adapters and filter out all reads that have a 31-mer match to [PhiX](https://emea.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html), which is commonly added to Illumina sequencing runs to monitor and/or improve overall run quality.
        
        ??? toggle "What are adapters and why do they need to be removed?"
            Adapters are manufactured oligonucleotide sequences attached to DNA fragments during the library preparation process. In Illumina sequencing, these adapter sequences are required for attaching reads to flow cells. You can read more about Illumina adapters [here](https://emea.support.illumina.com/bulletins/2020/06/illumina-adapter-portfolio.html). For genome analysis, it’s important to remove these sequences since they’re not actually from your sample. If you don’t remove them, the downstream analysis may be affected.

    ??? toggle "Read Quantification"

        There are two methods for read quantification to choose from: [`fastq-scan`](https://github.com/rpetit3/fastq-scan) (default) or [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Both quantify the forward and reverse reads in FASTQ files. In TheiaProk_Illumina_PE, they also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads. `fastqc` also provides a graphical visualization of the read quality.

    ??? toggle "Read Identification"
    
        Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

        Kraken2 is run on the set of raw reads, provided as input, as well as the set of clean reads that are resulted from the `read_QC_trim` workflow

        !!! info "Database-dependent"
            TheiaCoV automatically uses a viral-specific Kraken2 database.

        !!! techdetails "Kraken2 Technical Details"    
            
            |  | Links |
            | --- | --- |
            | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_kraken2.wdl) |
            | Software source code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
            | Software documentation | <https://github.com/DerrickWood/kraken2/wiki> |
            | Original publication | [Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |
       
    !!! techdetails "read_QC_trim Technical Details"
                
        | | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/workflows/wf_read_QC_trim.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_trimmomatic.wdl)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/quality_control/task_fastq_scan.wdl)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_midas.wdl)<br>[task_kraken2.wdl](https://github.com/theiagen/public_health_bacterial_genomics/blob/main/tasks/taxon_id/task_kraken2.wdl) |
        | Software source code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2)|
        | Software documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS); [Kraken2](https://github.com/DerrickWood/kraken2/wiki) |
        | Original publications | *[Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>*[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>*[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/)<br>*[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

??? task "`read_QC_trim_ONT`: Read Quality Trimming, Host Removal, and Identification _for ONT data_"

    `read_QC_trim` is a sub-workflow within TheiaCoV that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below.

    ??? toggle "Host removal"
    
        All reads of human origin **are removed**, including their mates, by using NCBI’s [**human read removal tool (HRRT)**](https://github.com/ncbi/sra-human-scrubber). 

        HRRT is based on the [SRA Taxonomy Analysis Tool](https://doi.org/10.1186/s13059-021-02490-0) and employs a k-mer database constructed of k-mers from Eukaryota derived from all human RefSeq records with any k-mers found in non-Eukaryota RefSeq records subtracted from the database.

        !!! techdetails "NCBI-Scrub Technical Details"
            
            |  | Links |
            | --- | --- |
            | Task | [task_ncbi_scrub.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/task_ncbi_scrub.wdl) |
            | Software source code | [NCBI Scrub on GitHub](https://github.com/ncbi/sra-human-scrubber) |
            | Software documentation | <https://github.com/ncbi/sra-human-scrubber/blob/master/README.md> |

    ??? toggle "Read quality filtering"
        
        Read filtering is performed using `artic guppyplex` which performs a quality check by filtering the reads by length to remove chimeric reads.
        
    ??? toggle "Read Identification"
    
        Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

        Kraken2 is run on the set of raw reads, provided as input, as well as the set of clean reads that are resulted from the `read_QC_trim` workflow

        !!! info "Database-dependent"
            TheiaCoV automatically uses a viral-specific Kraken2 database.

        !!! techdetails "Kraken2 Technical Details"    
            
            |  | Links |
            | --- | --- |
            | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_kraken2.wdl) |
            | Software source code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
            | Software documentation | <https://github.com/DerrickWood/kraken2/wiki> |
            | Original publication | [Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |
        
    !!! techdetails "read_QC_trim Technical Details"
        
        Each TheiaCoV workflow calls a sub-workflow listed below, which then calls the individual tasks:
        
        | Workflow | TheiaCoV_ONT |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim_ont.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim_ont.wdl) |
        | Tasks | [task_ncbi_scrub.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl#L68) (SE subtask)<br>[task_artic_guppyplex.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_artic_guppyplex.wdl)<br>[task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl#L3)|
        | Software source code | [NCBI Scrub on GitHub](https://github.com/ncbi/sra-human-scrubber)<br>[Artic on GitHub](https://github.com/artic-network/fieldbioinformatics)<br>[Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
        | Software documentation | [NCBI Scrub](<https://github.com/ncbi/sra-human-scrubber/blob/master/README.md>)<br>[Artic pipeline](https://artic.readthedocs.io/en/latest/?badge=latest)<br>[Kraken2](https://github.com/DerrickWood/kraken2/wiki) |
        | Original publications | [*STAT: a fast, scalable, MinHash-based *k*-mer tool to assess Sequence Read Archive next-generation sequence submissions](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02490-0)<br>*[Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0)  |

#### Assembly tasks

!!! tip ""
    Either one of these tasks is run depending on the organism and workflow type.

??? toggle "`ivar_consensus`: Alignment, Consensus, Variant Detection, and Assembly Statistics ==_for non-flu organisms in Illumina workflows_=="

    `ivar_consensus` is a sub-workflow within TheiaCoV that performs reference-based consensus assembly using the [iVar](https://andersen-lab.github.io/ivar/html/index.html) tool by Nathan Grubaugh from the Andersen lab.

    The following steps are performed as part of this sub-workflow:

    1. Cleaned reads are aligned to the appropriate reference genome (see also the [*organism-specific parameters and logic*](./theiacov.md#org-specific) section above) using [BWA](http://bio-bwa.sourceforge.net/) to generate a Binary Alignment Mapping (BAM) file.
    2. If `trim_primers` is set to true, primers will be removed using `ivar trim`.
        1.  General statistics about the remaining reads are calculated.
    3. The `ivar consensus` command is run to generate a consensus assembly.
    4. General statistics about the assembly are calculated..

    !!! techdetails "iVar Consensus Technical Details"    
        | Workflow | TheiaCoV_Illumina_PE & TheiaCoV_Illumina_SE |
        | --- | --- |
        | Sub-workflow | [wf_ivar_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_ivar_consensus.wdl) |
        | Tasks | [task_bwa.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/alignment/task_bwa.wdl)<br>[task_ivar_primer_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_ivar_primer_trim.wdl)<br>[task_assembly_metrics.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_assembly_metrics.wdl)<br>[task_ivar_variant_call.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/variant_detection/task_ivar_variant_call.wdl)<br>[task_ivar_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_ivar_consensus.wdl) |
        | Software source code | [BWA on GitHub](https://github.com/lh3/bwa), [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Software documentation | [BWA on SourceForge](https://bio-bwa.sourceforge.net/), [iVar on GitHub](https://andersen-lab.github.io/ivar/html/) |
        | Original publications | [*Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM](https://doi.org/10.48550/arXiv.1303.3997)<br>[*An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar](http://dx.doi.org/10.1186/s13059-018-1618-7) |

??? toggle "`artic_consensus`: Alignment, Primer Trimming, Variant Detection, and Consensus ==_for non-flu organisms in ONT & ClearLabs workflows_=="

    Briefly, input reads are aligned to the appropriate reference with [minimap2](https://github.com/lh3/minimap2) to generate a Binary Alignment Mapping ([BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map)) file. Primer sequences are then removed from the BAM file and a consensus assembly file is generated using the [Artic minion](https://artic.readthedocs.io/en/latest/commands/#basecaller) Medaka argument.

    !!! info ""
        Read-trimming is performed on raw read data generated on the ClearLabs instrument and thus not a required step in the TheiaCoV_ClearLabs workflow.

    General statistics about the assembly are generated with the `consensus_qc` task ([task_assembly_metrics.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_assembly_metrics.wdl)).

    !!! techdetails "Artic Consensus Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_artic_consensus.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_artic_consensus.wdl) |
        | Software source code | [Artic on GitHub](https://github.com/artic-network/fieldbioinformatics) |
        | Software documentation | [Artic pipeline](https://artic.readthedocs.io/en/latest/?badge=latest) |

??? toggle "`irma`: Assembly and Characterization ==_for flu in TheiaCoV_Illumina_PE & TheiaCoV_ONT_=="

    Cleaned reads are assembled using `irma` which does not use a reference due to the rapid evolution and high variability of influenza. `irma` also performs typing and subtyping as part of the assembly process.

    General statistics about the assembly are generated with the `consensus_qc` task ([task_assembly_metrics.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_assembly_metrics.wdl)).

    !!! techdetails "IRMA Technical Details" 
        |  | Links |
        | --- | --- |
        | Task | [task_irma.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_irma.wdl) |
        | Software documentation | [IRMA website](https://wonder.cdc.gov/amd/flu/irma/) |
        | **Original publications** | [*Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6) |

#### Organism-specific characterization tasks

The following table illustrates which characterization tools are run for the indicated organism.

|  | SARS-CoV-2 | MPXV | HIV | WNV | Influenza | RSV-A | RSV-B |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Pangolin | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Nextclade | ✅ | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ |
| VADR | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ | ✅ |
| Quasitools HyDRA | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ | ❌ |
| IRMA | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |
| Abricate | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |
| % Gene Coverage | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ❌ |
| Antiviral Detection | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |
| GenoFLU | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ | ❌ |

??? task "`pangolin`"

    Pangolin designates SARS-CoV-2 lineage assignments.
    
    !!! techdetails "Pangolin Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_pangolin.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/betacoronavirus/task_pangolin.wdl) |
        | Software source code | [Pangolin on GitHub](https://github.com/cov-lineages/pangolin) |
        | Software documentation | [Pangolin website](https://cov-lineages.org/resources/pangolin.html) |

??? task "`nextclade`"

    ["Nextclade is an open-source project for viral genome alignment, mutation calling, clade assignment, quality checks and phylogenetic placement."](https://docs.nextstrain.org/projects/nextclade/en/stable/)
    
    !!! techdetails "Nextclade Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_nextclade.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_nextclade.wdl#L63) |
        | Software source code | <https://github.com/nextstrain/nextclade> |
        | Software documentation | [Nextclade](https://docs.nextstrain.org/projects/nextclade/en/stable/) |
        | Original publication | [Nextclade: clade assignment, mutation calling and quality control for viral genomes.](https://doi.org/10.21105/joss.03773) |

??? task "`vadr`"

    VADR annotates and validates completed assembly files.

    !!! techdetails "VADR Technical Details"        
        
        |  | Links |
        | --- | --- |
        | Task | [task_vadr.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/advanced_metrics/task_vadr.wdl) |
        | Software source code | <https://github.com/ncbi/vadr> |
        | Software documentation | <https://github.com/ncbi/vadr/wiki> |
        | Original publication | For SARS-CoV-2: *[Faster SARS-CoV-2 sequence validation and annotation for GenBank using VADR](https://doi.org/10.1093/nargab/lqad002)*<br> For non-SARS_CoV-2: [*VADR: validation and annotation of virus sequence submissions to GenBank*](https://doi.org/10.1186/s12859-020-3537-3) |

??? task "`quasitools`"

    `quasitools` performs genome characterization for HIV.
    
    !!! techdetails "Quasitools Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_quasitools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/lentivirus/task_quasitools.wdl) |
        | Software source code | <https://github.com/phac-nml/quasitools/> |
        | Software documentation | [Quasitools HyDRA](https://phac-nml.github.io/quasitools/hydra/) |

??? task "`irma`"

    IRMA assigns types and subtype/lineages in addition to performing assembly of flu genomes. Please see the section above under "Assembly tasks" to find more information regarding this tool.
    
    !!! techdetails "IRMA Technical Details" 
        |  | Links |
        | --- | --- |
        | Task | [task_irma.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/assembly/task_irma.wdl) |
        | Software documentation | [IRMA website](https://wonder.cdc.gov/amd/flu/irma/) |
        | **Original publications** | [*Viral deep sequencing needs an adaptive approach: IRMA, the iterative refinement meta-assembler](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6) |

??? task "`abricate`"

    Abricate assigns types and subtype/lineages for flu samples
    
    !!! techdetails "Abricate Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_abricate.wdl (abricate_flu subtask)](https://github.com/theiagen/public_health_bioinformatics/blob/2dff853defc6ea540a058873f6fe6a78cc2350c7/tasks/gene_typing/drug_resistance/task_abricate.wdl#L59) |
        | Software source code | [ABRicate on GitHub](https://github.com/tseemann/abricate) |
        | Software documentation | [ABRicate on GitHub](https://github.com/tseemann/abricate) |

??? task "`gene_coverage`"

    This task calculates the percent of the gene covered above a minimum depth. By default, it runs for SARS-CoV-2 and MPXV, but if a bed file is provided with regions of interest, this task will be run for other organisms as well.

    !!! techdetails "Gene Coverage Technical Details"        
        |  | Links |
        | --- | --- |
        | Task | [task_gene_coverage.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_gene_coverage.wdl) |

??? task "`flu_antiviral_substitutions`"

    This sub-workflow determines which, if any, antiviral mutations are present in the sample. 
    
    The assembled HA, NA, PA, PB1 and PB2 segments are compared against [a list of known amino-acid substitutions associated with resistance](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/gene_typing/task_flu_antiviral_subs.wdl) to the antivirals  A_315675, compound_367, Favipiravir, Fludase, L_742_001, Laninamivir, Oseltamivir (tamiflu), Peramivir, Pimodivir, Xofluza, and Zanamivir. The list of known antiviral amino acid substitutions can be expanded via optional user input `antiviral_aa_subs` in the format "`NA:V95A,HA:I97V`", i.e. `Protein:AAPositionAA`. 

    !!! techdetails "Antiviral Substitutions Technical Details"        
        |  | Links |
        | --- | --- |
        | Workflow | [wf_influenza_antiviral_substitutions.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_influenza_antiviral_substitutions.wdl) |

??? task "`genoflu`"

    This sub-workflow determines the whole-genome genotype of an H5N1 flu sample.
    
    !!! techdetails "GenoFLU Technical Details"
        
        |  | Links |
        | --- | --- |
        | Task | [task_genoflu.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/species_typing/orthomyxoviridae/task_genoflu.wdl) |
        | Software source code | [GenoFLU on GitHub](https://github.com/USDA-VS/GenoFLU) |

### Outputs

[TheiaCoV Outputs (Not applicable to TheiaCoV_FASTA_Batch)](https://www.notion.so/4640fc5152f042bd91f74a0e308f7881?pvs=21)

!!! warning "Overwrite Warning"
    **TheiaCoV_FASTA_Batch_PHB** workflow will **output results to the set-level data table in addition to overwriting the Pangolin & Nextclade output columns in the sample-level data table**. Users can view the set-level workflow output TSV file called `"Datatable"` to view exactly which columns were overwritten in the sample-level data table.

[TheiaCoV_FASTA_Batch_PHB Outputs](https://www.notion.so/9ebe975169ec43c8a648b1121a187419?pvs=21)