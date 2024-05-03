version 1.0

import "../../tasks/gene_typing/plasmid_detection/task_tiptoft.wdl" as tiptoft_task
import "../../tasks/quality_control/read_filtering/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2
import "../../tasks/utilities/task_kmc.wdl" as kmc_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads with (1) fastq_scan, (2) nanoplot, (3) kmc, (4) rasusa downsampling, (5) tiptoft plasmid detection, and (6) nanoq filtering"
  }
  input {
    String samplename
    File read1
    Int? genome_length

    String? workflow_series

    # theiacov inputs
    Int? min_length
    Int? max_length
    String? run_prefix

    # kraken inputs
    String? target_organism

    # rasusa downsampling
    Float downsampling_coverage = 150
  }
  if ("~{workflow_series}" == "theiacov") {
    call ncbi_scrub.ncbi_scrub_se {
      input:
        read1 = read1,
        samplename = samplename
    }
    call artic_guppyplex.read_filtering {
      input:
        read1 = ncbi_scrub_se.read1_dehosted,
        samplename = samplename,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix
    }
    call kraken2.kraken2_theiacov as kraken2_raw {
      input:
        samplename = samplename,
        read1 = read1,
        target_organism = target_organism
    }
    call kraken2.kraken2_parse_classified as kraken2_recalculate_abundances_raw {
      input:
        samplename = samplename,
        kraken2_report = kraken2_raw.kraken_report,
        kraken2_classified_report = kraken2_raw.kraken2_classified_report,
        target_organism = target_organism
    }  
    call kraken2.kraken2_theiacov as kraken2_dehosted {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_se.read1_dehosted,
        target_organism = target_organism
    }
    call kraken2.kraken2_parse_classified as kraken2_recalculate_abundances_dehosted {
      input:
        samplename = samplename,
        kraken2_report = kraken2_dehosted.kraken_report,
        kraken2_classified_report = kraken2_dehosted.kraken2_classified_report,
        target_organism = target_organism
    } 
  }
  if ("~{workflow_series}" == "theiaprok") {
    # kmc for genome size estimation
    call kmc_task.kmc {
      input:
        read1 = read1,
        samplename = samplename
    }

    Int kmc_est_genome_length = if kmc.est_genome_length > 10000 then 10000 else kmc.est_genome_length

    # rasusa for random downsampling
    call rasusa_task.rasusa {
      input:
        read1 = read1,
        samplename = samplename,
        coverage = downsampling_coverage,
        genome_length = select_first([genome_length, kmc_est_genome_length])
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
    String? kraken_target_organism_name = kraken2_raw.kraken_target_organism
    
    # kraken outputs raw
    Float? kraken_human = kraken2_recalculate_abundances_raw.percent_human
    Float? kraken_sc2 = kraken2_recalculate_abundances_raw.percent_sc2
    String? kraken_target_organism = kraken2_recalculate_abundances_raw.percent_target_organism
    File? kraken_report = kraken2_recalculate_abundances_raw.kraken_report
    
    # kraken outputs dehosted
    Float? kraken_human_dehosted = kraken2_recalculate_abundances_dehosted.percent_human
    Float? kraken_sc2_dehosted = kraken2_recalculate_abundances_dehosted.percent_sc2
    String? kraken_target_organism_dehosted = kraken2_recalculate_abundances_dehosted.percent_target_organism
    File? kraken_report_dehosted = kraken2_recalculate_abundances_dehosted.kraken_report
   
    # theiaprok outputs
    # kmc outputs
    Int? est_genome_length = kmc_est_genome_length
    File? kmc_kmer_stats = kmc.kmer_stats
    String? kmc_version = kmc.kmc_version
    
    # nanoq outputs
    File read1_clean = select_first([nanoq.filtered_read1, read_filtering.read1_clean])
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