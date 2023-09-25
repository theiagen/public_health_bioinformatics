version 1.0

import "../../tasks/quality_control/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/quality_control/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/task_kraken2.wdl" as kraken2
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task
import "../../tasks/utilities/task_kmc.wdl" as kmc_task
import "../../tasks/gene_typing/task_tiptoft.wdl" as tiptoft_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads with (1) fastq_scan, (2) nanoplot, (3) kmc, (4) rasusa downsampling, (5) tiptoft plasmid detection, and (6) nanoq filtering"
  }
  input {
    String samplename
    File read1
    Int? genome_size

    String? workflow_series

    # theiacov inputs
    Int? min_length
    Int? max_length
    String? run_prefix

    # kraken inputs
    String? target_org
  }
  if ("~{workflow_series}" == "theiacov") {
    call ncbi_scrub.ncbi_scrub_se {
      input:
        read1 = read1,
        samplename = samplename
    }
    call artic_guppyplex.read_filtering {
      input:
        demultiplexed_reads = ncbi_scrub_se.read1_dehosted,
        samplename = samplename,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix
    }
    call kraken2.kraken2_theiacov as kraken2_raw {
      input:
        samplename = samplename,
        read1 = read1,
        target_org = target_org
    }  
    call kraken2.kraken2_theiacov as kraken2_dehosted {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_se.read1_dehosted,
        target_org = target_org
    }
  }
  if ("~{workflow_series}" == "theiaprok") {
    # kmc for genome size estimation
    call kmc_task.kmc {
      input:
        read1 = read1,
        samplename = samplename
    }
    # rasusa for random downsampling
    call rasusa_task.rasusa {
      input:
        read1 = read1,
        samplename = samplename,
        coverage = 150,
        genome_size = select_first([genome_size, kmc.est_genome_size])
    }
    # tiptoft for plasmid detection
    call tiptoft_task.tiptoft {
      input:
        read1 = read1,
        samplename = samplename
    }  
    # nanoq/filtlong (default min length 500)
    call nanoq_task.nanoq {
      input:
        read1 = rasusa.read1_subsampled,
        samplename = samplename
    }
  }
  output { 
    # theiacov outputs
    # ncbi scrub outputs
    File? read1_dehosted = ncbi_scrub_se.read1_dehosted
    
    # kraken outputs
    String? kraken_version = kraken2_raw.version
    String? kraken_target_org_name = kraken2_raw.kraken_target_org
    
    # kraken outputs raw
    Float? kraken_human = kraken2_raw.percent_human
    Float? kraken_sc2 = kraken2_raw.percent_sc2
    String? kraken_target_org = kraken2_raw.percent_target_org
    File? kraken_report = kraken2_raw.kraken_report
    
    # kraken outputs dehosted
    Float? kraken_human_dehosted = kraken2_dehosted.percent_human
    Float? kraken_sc2_dehosted = kraken2_dehosted.percent_sc2
    String? kraken_target_org_dehosted = kraken2_dehosted.percent_target_org
    File? kraken_report_dehosted = kraken2_dehosted.kraken_report
   
    # theiaprok outputs
    # kmc outputs
    Int? est_genome_size = kmc.est_genome_size
    File? kmc_kmer_stats = kmc.kmer_stats
    String? kmc_version = kmc.kmc_version
    
    # nanoq outputs
    File read1_clean = select_first([nanoq.filtered_read1, read_filtering.filtered_reads])
    String? nanoq_version = nanoq.version

    # rasusa outputs
    String? rasusa_version = rasusa.rasusa_version  

    # tiptoft outputs
    File? tiptoft_plasmid_replicon_fastq = tiptoft.tiptoft_plasmid_replicon_fastq
    File? tiptoft_result_tsv = tiptoft.tiptoft_tsv
    String? tiptoft_plasmid_replicon_genes = tiptoft.plasmid_replicon_genes
    String? tiptoft_version = tiptoft.tiptoft_version
  }
}