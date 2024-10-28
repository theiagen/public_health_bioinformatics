# TheiaMeta Panel

## Quick Facts

| **Workflow Type** | **Applicable Kingdom** | **Last Known Changes** | **Command-line Compatibility** | **Workflow Level** |
|---|---|---|---|---|
| [Genomice Characterization](../../workflows_overview/workflows_type.md/#genomic_characterization) | [Viral](../../workflows_overview/workflows_kingdom.md/#viral) | PHB v2.X.X | Yes | Sample-level |

## TheiaMeta_Panel_Illumina_PE_PHB

TheiaMeta_Panel_Illumina_PE was created initially for the [Illumina Viral Surveillance Panel](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/viral-surveillance-panel.html)[^1]; however, it can be used for any panel that is sequenced using Illumina paired-end reads if the appropriate taxon IDs are provided. TheiaMeta_Panel performs taxonomic binning, and then assembles the bins into contigs. If the contigs are associated with a supported organism, genomic characterization will be performed.

[^1]: We are not affiliated with Illumina, Inc. The mention of the Illumina Viral Surveillance Panel is for informational purposes only.

??? toggle "**What organisms and taxon IDs are identified by default?**"
    The Illumina VSP panel contains over 224 viral species, of which 163 can be identified in the default Kraken2 viral database.

    Accordingly, the following 163 taxon IDs are used by default in TheiaMeta_Panel_Illumina_PE. Feel free to search this table to see if your organism of interest is included.

    <div class="searchable-table" markdown="1">

    | **Taxon ID** | **Organism Name in Illumina VSP Panel** |
    |---|---|
    | 10804 | Adeno-associated virus 2 (AAV2) |
    | 1313215 | Aichi virus 1 (AiV-A1) |
    | 2849717  | Aigai virus (AIGV) |
    | 1980456 | Andes virus (ANDV) |
    | 1424613 | Anjozorobe virus (ANJV) |
    | 90961 | Australian bat lyssavirus (ABLV) |
    | 3052470 | Bayou virus (BAYV) |
    | 3052490 | Black Creek Canal virus (BCCV) |
    | 2010960 | Bombali virus (BOMV) |
    | 1618189 | Bourbon virus (BRBV) |
    | 565995 | Bundibugyo virus (BDBV) |
    | 80935 | Cache Valley virus (CVV) |
    | 35305 | California encephalitis virus (CEV) |
    | 1221391 | Cedar virus (CedV) |
    | 3052302 | Chapare virus (CHAPV) |
    | 37124 | Chikungunya virus (CHIKV) |
    | 169173 | Choclo virus (CHOV) |
    | 46839 | Colorado tick fever virus (CTFV) |
    | 138948 | Coxsackievirus A |
    | 138949 | Coxsackievirus B |
    | 3052518 | Crimean-Congo hemorrhagic fever virus (CCHFV) |
    | 11053 | Dengue Virus 1 |
    | 11060 | Dengue Virus 2 |
    | 11069 | Dengue Virus 3 |
    | 11070 | Dengue Virus 4 |
    | 3052477 | Dobrava virus (DOBV) |
    | 38767 | Duvenhage virus (DUVV) |
    | 11021 | Eastern equine encephalitis virus (EEEV) |
    | 138951 | Enterovirus D |
    | 10376 | Epstein-Barr virus (EBV) |
    | 57482 | European bat lyssavirus 1 |
    | 57483 | European bat lyssavirus 2 |
    | 2847089 | Ghana virus (GhV) |
    | 3052307 | Guanarito virus (GTOV) |
    | 3052480 | Hantaan virus (HTNV) |
    | 1216928 | Heartland virus (HRTV) |
    | 3052223 | Hendra virus (HeV) |
    | 12092 | Hepatitis A virus (HAV) |
    | 3052230 | Hepatitis C virus (HCV) |
    | 12475 | Hepatitis D virus (HDV) |
    | 10298 | Herpes simplex virus 1 (HSV1) |
    | 129875 | Human adenovirus A |
    | 108098 | Human adenovirus B |
    | 129951 | Human adenovirus C |
    | 130310 | Human adenovirus D |
    | 130308 | Human adenovirus E |
    | 130309 | Human adenovirus F |
    | 536079 | Human adenovirus G |
    | 11137 | Human coronavirus 229E (HCoV_229E) |
    | 290028 | Human coronavirus HKU1 (HCoV_HKU1) |
    | 277944 | Human coronavirus NL63 (HCoV_NL63) |
    | 31631 | Human coronavirus OC43 (HCoV_OC43) |
    | 10359 | Human cytomegalovirus (HCMV) |
    | 11676 | Human immunodeficiency virus 1 (HIV-1) |
    | 11709 | Human immunodeficiency virus 2 (HIV-2) |
    | 162145 | Human metapneumovirus (HMPV) |
    | 333760 | Human papillomavirus 16 (HPV16; high-risk) |
    | 333761 | Human papillomavirus 18 (HPV18; high-risk) |
    | 333762 | Human papillomavirus 26 (HPV26) |
    | 12730 | Human parainfluenza virus 1 (HPIV-1) |
    | 2560525 | Human parainfluenza virus 2 (HPIV-2) |
    | 11216 | Human parainfluenza virus 3 (HPIV-3) |
    | 2560526  | Human parainfluenza virus 4 (HPIV-4) |
    | 1803956  | Human parechovirus (HPeV) |
    | 10798  | Human parvovirus B19 (B19V) |
    | 746830 | Human polyomavirus 6 (HPyV6) |
    | 746831 | Human polyomavirus 7 (HPyV7) |
    | 943908 | Human polyomavirus 9 (HPyV9) |
    | 208893 | Human respiratory syncytial virus A (HRSV-A) |
    | 114727 | Influenza A virus (H1N1) |
    | 114729 | Influenza A virus (H2N2) |
    | 119210 | Influenza A virus (H3N2) |
    | 102793 | Influenza A virus (H5N1) |
    | 333278 | Influenza A virus (H7N9) |
    | 102796 | Influenza A virus (H9N2) |
    | 11520 | Influenza B virus |
    | 11552 | Influenza C virus |
    | 35511 | Jamestown Canyon virus (JCV) |
    | 11072 | Japanese encephalitis virus (JEV) |
    | 10632 | JC polyomavirus (JCPyV) |
    | 2169991 | Junin virus (JUNV) |
    | 1891764 | KI polyomavirus (KIPyV) |
    | 33743 | Kyasanur Forest disease virus (KFDV) |
    | 11577 | La Crosse virus (LACV) |
    | 38766 | Lagos bat virus (LBV) |
    | 3052489 | Laguna Negra virus (LANV) |
    | 3052310 | Lassa virus (LASV) |
    | 1965344 | LI polyomavirus (LIPyV) |
    | 3052148 | Lloviu virus (LLOV) |
    | 3052314 | Lujo virus (LUJV) |
    | 3052303 | Lymphocytic choriomeningitis virus (LCMV) |
    | 3052317 | Machupo virus (MACV) |
    | 1239565 | Mamastrovirus 1 (MAstV1) |
    | 1239570 | Mamastrovirus 6 (MAstV6) |
    | 1239573 | Mamastrovirus 9 (MAstV9) |
    | 238817 | Maporal virus (MAPV) |
    | 3052505 | Marburg virus (MARV) |
    | 59301 | Mayaro virus (MAYV) |
    | 11234 | Measles virus (MV) |
    | 152219 | Menangle virus (MenV) |
    | 493803 | Merkel cell polyomavirus (MCPyV) |
    | 1335626 | Middle East respiratory syndrome-related coronavirus (MERS-CoV) |
    | 1474807 | Mojiang virus (MojV) |
    | 12538 | Mokola virus (MOKV) |
    | 10244 | Monkeypox virus (MPV) |
    | 2560602 | Mumps virus (MuV) |
    | 11079 | Murray Valley encephalitis virus (MVEV) |
    | 1203539 | MW polyomavirus (MWPyV) |
    | 1497391 | New Jersey polyomavirus (NJPyV) |
    | 3052225 | Nipah virus (NiV) |
    | 142786 | Norovirus |
    | 12542 | Omsk hemorrhagic fever virus (OHFV) |
    | 2169701 | Onyong-nyong virus (ONNV) |
    | 118655 | Oropouche virus (OROV) |
    | 138950 | Poliovirus |
    | 11083 | Powassan virus (POWV) |
    | 11587 | Punta Toro virus (PTV) |
    | 3052493 | Puumala virus (PUUV) |
    | 11292 | Rabies virus (RABV) |
    | 186539 | Reston virus (RESTV) |
    | 147711 | Rhinovirus A (RV-A) |
    | 147712 | Rhinovirus B (RV-B) |
    | 463676 | Rhinovirus C (RV-C) |
    | 11588 | Rift Valley fever virus (RVFV) |
    | 11029 | Ross River virus (RRV) |
    | 28875 | Rotavirus A (RVA) |
    | 28876 | Rotavirus B (RVB) |
    | 36427 | Rotavirus C (RVC) |
    | 1348384 | Rotavirus H (RVH) |
    | 11041 | Rubella virus (RuV) |
    | 2907957 | Sabia virus (SBAV) |
    | 1330524 | Salivirus A (SaV-A) |
    | 3052496 | Sangassou virus (SANGV) |
    | 95341 | Sapovirus |
    | 11033 | Semliki Forest virus (SFV) |
    | 3052498 | Seoul virus (SEOV) |
    | 2901879 | Severe acute respiratory syndrome coronavirus (SARS-CoV) |
    | 2697049 | Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) |
    | 1003835 | Severe fever with thrombocytopenia syndrome virus (SFTSV) |
    | 1891767 | Simian virus 40 (SV40) |
    | 3052499 | Sin nombre virus (SNV) |
    | 11034 | Sindbis virus (SINV) |
    | 11580 | Snowshoe hare virus (SSHV) |
    | 1452514 | Sosuga virus (SoRV) |
    | 11080 | St. Louis encephalitis virus (SLEV) |
    | 1277649 | STL polyomavirus (STLPyV) |
    | 186540 | Sudan virus (SUDV) |
    | 1608084 | Tacheng tick virus 2 (TcTV-2) |
    | 45270 | Tahyna virus (TAHV) |
    | 186541 | Tai Forest virus (TAFV) |
    | 11084 | Tick-borne encephalitis virus (TBEV) |
    | 68887 | Torque teno virus (TTV) |
    | 862909 | Trichodysplasia spinulosa-associated polyomavirus (TSPyV) |
    | 3052503 | Tula virus (TULV) |
    | 64286 | Usutu virus (USUV) |
    | 10255 | Variola virus (VARV) |
    | 11036 | Venezuelan equine encephalitis virus (VEEV) |
    | 11082 | West Nile virus (WNV) |
    | 11039 | Western equine encephalitis virus (WEEV) |
    | 440266 | WU polyomavirus (WUPyV) |
    | 11089 | Yellow fever virus (YFV) |
    | 186538 | Zaire ebolavirus(EBOV) |
    | 64320 | Zika virus (ZIKV) |

    </div>

!!! tip "Make your own list of taxon IDs"
    You may want to make your own list of taxon IDs if you know your sample is likely to contain a specific organism or group of organisms. You can find taxon IDs in the [NCBI Taxonomy Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).

    In Terra, you provide your created list of taxon IDs as an array of integers for the `taxon_ids` optional input variable, like this: `[1, 2, 3, 4, 5]`. Just replace the numbers in this example with the taxon IDs you want to use.

### Inputs

<div class="searchable-table" markdown="1">

| **Terra Task Name** | **Variable** | **Type** | **Description** | **Default Value** | **Terra Status** |
|---|---|---|---|---|---|
| theiameta_panel_illumina_pe | **read1** | File | The forward Illumina read in FASTQ file format (compression optional)  | | Required |
| theiameta_panel_illumina_pe | **read2** | File | The reverse Illumina read in FASTQ file format (compression optional) | | Required |
| theiameta_panel_illumina_pe | **samplename** | String | The name of the sample being analyzed | | Required |
| fastq_scan_binned | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| fastq_scan_binned | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| fastq_scan_binned | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1 | Optional |
| fastq_scan_binned | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| gather_scatter | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| gather_scatter | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| gather_scatter | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16 | Optional |
| gather_scatter | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| kraken2 | **classified_out** | String | Allows user to rename the classified FASTQ files output. Must include .fastq as the suffix | classified#.fastq | Optional |
| kraken2 | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| kraken2 | **disk_size** | Int | GB of storage to request for VM used to run the kraken2 task. Increase this when using large (>30GB kraken2 databases such as the "k2_standard" database) | 100 | Optional |
| kraken2 | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/kraken2:2.1.2-no-db | Optional |
| kraken2 | **kraken2_args** | String | Allows a user to supply additional kraken2 command-line arguments |  | Optional |
| kraken2 | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| kraken2 | **unclassified_out** | String | Allows user to rename unclassified FASTQ files output. Must include .fastq as the suffix | unclassified#.fastq | Optional |
| krakentools | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| krakentools | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| krakentools | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/krakentools:d4a2fbe| Optional |
| krakentools | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| metaspades | **cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| metaspades | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| metaspades | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/theiagen/metaspades:3.15.3 | Optional |
| metaspades | **kmers** | String | The k-mer list to use; if not provided, the value is automatically set | | Optional |
| metaspades | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 32 | Optional |
| metaspades | **metaspades_opts** | String | Additional arguments to pass on to the metaspades command | | Optional |
| metaspades | **phred_offset** | Int | The PHRED quality offset of the input reads; can be either 33 or 64 | 33 | Optional |
| minimap2_assembly_correction | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| minimap2_assembly_correction | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| minimap2_assembly_correction | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/minimap2:2.22 | Optional |
| minimap2_assembly_correction | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| morgana_magic | **abricate_flu_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| morgana_magic | **abricate_flu_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| morgana_magic | **abricate_flu_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/abricate:1.0.1-insaflu-220727 | Optional |
| morgana_magic | **abricate_flu_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| morgana_magic | **abricate_flu_mincov** | Int | Minimum DNA % coverage | 60 | Optional |
| morgana_magic | **abricate_flu_minid** | Int | Minimum DNA % identity | 70 | Optional |
| morgana_magic | **assembly_metrics_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| morgana_magic | **assembly_metrics_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| morgana_magic | **assembly_metrics_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15 | Optional |
| morgana_magic | **assembly_metrics_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| morgana_magic | **cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| morgana_magic | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| morgana_magic | **docker** | String | The Docker container to use for the task | ngolin | Optional |
| morgana_magic | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| morgana_magic | **genoflu_cpu** | Int | Number of CPUs to allocate to the task | 1 | Optional |
| morgana_magic | **genoflu_cross_reference** | File | An Excel file to cross-reference BLAST findings; probably useful if novel genotypes are not in the default file used by genoflu.py | | Optional |
| morgana_magic | **genoflu_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 25 | Optional |
| morgana_magic | **genoflu_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/genoflu:1.03 | Optional |
| morgana_magic | **genoflu_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 2 | Optional |
| morgana_magic | **irma_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| morgana_magic | **irma_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| morgana_magic | **irma_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/cdcgov/irma:v1.1.5 | Optional |
| morgana_magic | **irma_keep_ref_deletions** | Boolean | True/False variable that determines if sites missed during read gathering should be deleted by ambiguation. | TRUE | Optional |
| morgana_magic | **irma_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 16 | Optional |
| morgana_magic | **nextclade_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| morgana_magic | **nextclade_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| morgana_magic | **nextclade_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/nextstrain/nextclade:3.3.1 | Optional |
| morgana_magic | **nextclade_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| morgana_magic | **nextclade_output_parser_cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| morgana_magic | **nextclade_output_parser_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 50 | Optional |
| morgana_magic | **nextclade_output_parser_docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/python/python:3.8.18-slim | Optional |
| morgana_magic | **nextclade_output_parser_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 4 | Optional |
| morgana_magic | **pangolin_cpu** | Int | Number of CPUs to allocate to the task | 4 | Optional |
| morgana_magic | **pangolin_disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| morgana_magic |  **pangolin_docker_image** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pangolin:4.3.1-pdata-1.29 | Optional |
| morgana_magic | **pangolin_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| pilon | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| pilon | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| pilon | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/pilon:1.24 | Optional |
| pilon | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| quast | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| quast | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| quast | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2 | Optional |
| quast | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| read_QC_trim | **adapters** | File | A file containing the sequence of the adapters used during library preparation, used in the BBDuk task |  | Optional |
| read_QC_trim | **bbduk_memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| read_QC_trim | **call_kraken** | Boolean | Set to true to launch Kraken2; if true, you must provide a kraken_db | FALSE | Optional |
| read_QC_trim | **call_midas** | Boolean | Set to true to launch Midas | TRUE | Optional |
| read_QC_trim | **fastp_args** | String | Additional arguments to pass to fastp | "--detect_adapter_for_pe -g -5 20 -3 20 | Optional |
| read_QC_trim | **midas_db** | File | Midas database file | gs://theiagen-large-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz | Optional |
| read_QC_trim | **phix** | File | A file containing the phix used during Illumina sequencing; used in the BBDuk task |  | Optional |
| read_QC_trim | **read_processing** | String | Read trimming software to use, either "trimmomatic" or "fastp" | trimmomatic | Optional |
| read_QC_trim | **read_qc** | String | Allows the user to decide between fastq_scan (default) and fastqc for the evaluation of read quality. | fastq_scan | Optional |
| read_QC_trim | **trim_min_length** | Int | The minimum length of each read after trimming | 75 | Optional | 
| read_QC_trim | **trim_primers** | Boolean | A True/False option that determines if primers should be trimmed. | TRUE | Optional |
| read_QC_trim | **trim_quality_min_score** | Int | The minimum quality score to keep during trimming | 30 | Optional |
| read_QC_trim | **trim_window_size** | Int | Specifies window size for trimming (the number of bases to average the quality across) | 4 | Optional |
| read_QC_trim | **trimmomatic_args** | String | Additional arguments to pass to trimmomatic | -phred33 | Optional |
| sort_bam_assembly_correction | **cpu** | Int | Number of CPUs to allocate to the task | 2 | Optional |
| sort_bam_assembly_correction | **disk_size** | Int | Amount of storage (in GB) to allocate to the task | 100 | Optional |
| sort_bam_assembly_correction | **docker** | String | The Docker container to use for the task | us-docker.pkg.dev/general-theiagen/staphb/samtools:1.17 | Optional |
| sort_bam_assembly_correction | **memory** | Int | Amount of memory/RAM (in GB) to allocate to the task | 8 | Optional |
| theiameta_panel_illumina_pe | **kraken2_db** | File | A Kraken2 database in .tar.gz format | gs://theiagen-large-public-files-rp/terra/databases/kraken2/k2_viral_20240112.tar.gz | Optional |
| theiameta_panel_illumina_pe | **minimum_read_number** | Int | The minimum number of reads in order to attempt assembly on a bin of reads | 1000 | Optional |
| theiameta_panel_illumina_pe | **taxon_ids** | Array[Int] | The taxon IDs to be used for taxonomic binning. By default, this array uses the taxon IDs listed above that are intended for the Illumina VSP panel | Illumina VSP panel (see above toggle) | Optional |
| version_capture | **docker** | String | The Docker container to use for the task | "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0" | Optional |
| version_capture | **timezone** | String | Set the time zone to get an accurate date of analysis (uses UTC by default) | | Optional |

</div>

### Workflow Tasks

??? task "`read_QC_trim`: Read Quality Trimming, Adapter Removal, Quantification, and Identification"

    `read_QC_trim` is a sub-workflow within TheiaMeta that removes low-quality reads, low-quality regions of reads, and sequencing adapters to improve data quality. It uses a number of tasks, described below.

    **Read quality trimming**

    Either `trimmomatic` or `fastp` can be used for read-quality trimming. Trimmomatic is used by default. Both tools trim low-quality regions of reads with a sliding window (with a window size of `trim_window_size`), cutting once the average quality within the window falls below `trim_quality_trim_score`. They will both discard the read if it is trimmed below `trim_min_length`. 

    By default, the trim_min_length is set to 75 bp. This is likely _too high_ for data generated using the Illumina VSP panel. We recommend setting this parameter to `50` in this case.

    If fastp is selected for analysis, fastp also implements the additional read-trimming parameters indicated below:

    | **Parameter** | **Explanation** |
    | --- | --- |
    | -g | enables polyG tail trimming |
    | -5 20 | enables read end-trimming |
    | -3 20 | enables read end-trimming |
    | --detect_adapter_for_pe | enables adapter-trimming **only for paired-end reads** |

    **Adapter removal**

    The `BBDuk` task removes adapters from sequence reads. To do this:

    - [Repair](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/repair-guide/) from the [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) package reorders reads in paired fastq files to ensure the forward and reverse reads of a pair are in the same position in the two fastq files.
    - [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)  (*"Bestus Bioinformaticus" Decontamination Using Kmers*) is then used to trim the adapters and filter out all reads that have a 31-mer match to [PhiX](https://emea.illumina.com/products/by-type/sequencing-kits/cluster-gen-sequencing-reagents/phix-control-v3.html), which is commonly added to Illumina sequencing runs to monitor and/or improve overall run quality.
    
    ??? toggle "What are adapters and why do they need to be removed?"
        Adapters are manufactured oligonucleotide sequences attached to DNA fragments during the library preparation process. In Illumina sequencing, these adapter sequences are required for attaching reads to flow cells. You can read more about Illumina adapters [here](https://emea.support.illumina.com/bulletins/2020/06/illumina-adapter-portfolio.html). For genome analysis, it's important to remove these sequences since they're not actually from your sample. If you don't remove them, the downstream analysis may be affected.
        
    **Read Quantification**

    There are two methods for read quantification to choose from: [`fastq-scan`](https://github.com/rpetit3/fastq-scan) (default) or [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Both quantify the forward and reverse reads in FASTQ files. In paired-end workflows, they also provide the total number of read pairs. This task is run once with raw reads as input and once with clean reads as input. If QC has been performed correctly, you should expect **fewer** clean reads than raw reads. `fastqc` also provides a graphical visualization of the read quality in an HTML file.

    **Read Identification (optional)**

    The `MIDAS` task is for the identification of reads to detect contamination with non-target taxa. This task is optional and turned off by default. It can be used by setting the `call_midas` input variable to `true`.

    The MIDAS reference database, located at **`gs://theiagen-large-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz`**, is provided as the default. It is possible to provide a custom database. More information is available [here](https://github.com/snayfach/MIDAS/blob/master/docs/ref_db.md).

    ??? toggle "How are the MIDAS output columns determined?"
        
        Example MIDAS report in the ****`midas_report` column:
        
        | species_id | count_reads | coverage | relative_abundance |
        | --- | --- | --- | --- |
        | Salmonella_enterica_58156 | 3309 | 89.88006645 | 0.855888033 |
        | Salmonella_enterica_58266 | 501 | 11.60606061 | 0.110519371 |
        | Salmonella_enterica_53987 | 99 | 2.232896237 | 0.021262881 |
        | Citrobacter_youngae_61659 | 46 | 0.995216227 | 0.009477003 |
        | Escherichia_coli_58110 | 5 | 0.123668877 | 0.001177644 |
        
        MIDAS report column descriptions:
        
        - species_id: species identifier
        - count_reads: number of reads mapped to marker genes
        - coverage: estimated genome-coverage (i.e. read-depth) of species in metagenome
        - relative_abundance: estimated relative abundance of species in metagenome
  
    !!! techdetails "read_QC_trim Technical Details"
                
        |  | Links |
        | --- | --- |
        | Sub-workflow | [wf_read_QC_trim.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/workflows/utilities/wf_read_QC_trim.wdl) |
        | Tasks | [task_fastp.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_fastp.wdl)<br>[task_trimmomatic.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_trimmomatic.wdl)<br>[task_bbduk.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/read_filtering/task_bbduk.wdl)<br>[task_fastq_scan.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/quality_control/basic_statistics/task_fastq_scan.wdl)<br>[task_midas.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_midas.wdl) |
        | Software Source Code | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](https://github.com/usadellab/Trimmomatic); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS)|
        | Software Documentation | [fastp](https://github.com/OpenGene/fastp); [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic); [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/); [fastq-scan](https://github.com/rpetit3/fastq-scan); [MIDAS](https://github.com/snayfach/MIDAS) |
        | Original Publication(s) | [Trimmomatic: a flexible trimmer for Illumina sequence data](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4103590/)<br>[fastp: an ultra-fast all-in-one FASTQ preprocessor](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234?login=false)<br>[An integrated metagenomics pipeline for strain profiling reveals novel patterns of bacterial transmission and biogeography](https://pubmed.ncbi.nlm.nih.gov/27803195/) |

??? task "`kraken2`: Taxonomic Classification"

    Kraken2 is a bioinformatics tool originally designed for metagenomic applications. It has additionally proven valuable for validating taxonomic assignments and checking contamination of single-species (e.g. bacterial isolate, eukaryotic isolate, viral isolate, etc.) whole genome sequence data.

    Kraken2 is run on the clean reads that result from the `read_QC_trim` subworkflow. By default, the Kraken2 database is set to the `k2_viral_20240112` database, located at `"gs://theiagen-large-public-files-rp/terra/databases/kraken2/k2_viral_20240112.tar.gz"`.

    !!! info "Database-dependent"
        The Kraken2 software is database-dependent and **taxonomic assignments are highly sensitive to the database used**. An appropriate database should contain the expected organism(s) (e.g. _Escherichia coli_) and other taxa that may be present in the reads (e.g. _Citrobacter freundii_, a common contaminant).

    !!! techdetails "Kraken2 Technical Details"    
        |  | Links |
        | --- | --- |
        | Task | [task_kraken2.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/contamination/task_kraken2.wdl) |
        | Software Source Code | [Kraken2 on GitHub](https://github.com/DerrickWood/kraken2/) |
        | Software Documentation | <https://github.com/DerrickWood/kraken2/wiki> |
        | Original Publication(s) | [Improved metagenomic analysis with Kraken 2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1891-0) |

??? task "`extract_kraken_reads` from KrakenTools: Read Binning"
    KrakenTools is a collection of scripts that can be used to help downstream analysis of Kraken2 results. In particular, this task uses the `extract_kraken_reads` script, which extracts reads classified at any user-specified taxonomy IDs. All parent and children reads of the specified taxonomic ID are also extracted.

    !!! techdetails "KrakenTools Technical Details"
        |  | Links |
        | --- | --- |
        | Task | [task_kraken_tools.wdl](https://github.com/theiagen/public_health_bioinformatics/blob/main/tasks/taxon_id/task_krakentools.wdl)
        | Software Source Code | [KrakenTools on GitHub](https://github.com/jenniferlu717/KrakenTools) |
        | Software Documentation | [KrakenTools on GitHub](https://github.com/jenniferlu717/KrakenTools) |
        | Original Publication | [Metagenome analysis using the Kraken software suite](https://doi.org/10.1038/s41596-022-00738-y) |

### Outputs

<div class="searchable-table" markdown="1">

| **Variable** | **Type** | **Description** |
|---|---|---|
| identified_organisms | Array[String] | A list of organisms that were able to be identified in the sample with the specified Kraken2 database |
| kraken2_classified_report | File | Standard Kraken2 output report. TXT filetype, but can be opened in Excel as a TSV file |
| kraken2_database | String | The name of the database used to run Kraken2 |
| kraken2_docker | String | Docker image used to run kraken2 |
| kraken2_report | File | Text document describing taxonomic prediction of every FASTQ record. This file can be very large and cumbersome to open and view |
| kraken2_version | String | The version of Kraken2 used in the analysis |
| results_by_taxon_tsv | File | A TSV file that contains the results for every taxon ID provided in the taxon_ids input variable that had reads identified; characterization (if applicable) and basic statistics regarding read count, assembly generation (if applicable), and general quality, are also associated with each bin |
| theiameta_panel_illumina_pe_analysis_date | String | Date the workflow was run |
| theiameta_panel_illumina_pe_version | String | Version of PHB used to run the workflow |

</div>

#### The `results_by_taxon_tsv` Output File

This file contains the 

## References

