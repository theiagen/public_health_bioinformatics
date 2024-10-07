version 1.0

import "../../tasks/gene_typing/plasmid_detection/task_tiptoft.wdl" as tiptoft_task
import "../../tasks/quality_control/read_filtering/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads with (1) fastq_scan, (2) nanoplot, (3) rasusa downsampling, (4) tiptoft plasmid detection, and (5) nanoq filtering"
  }
  input {
    String samplename
    File read1

    # kmc has been observed to be unreliable for genome length estimation, so we are now using a fixed value
    # setting this to be 5Mb which is around .7Mb greater than the mean genome length of bacteria (based on https://github.com/CDCgov/phoenix/blob/717d19c19338373fc0f89eba30757fe5cfb3e18a/assets/databases/NCBI_Assembly_stats_20240124.txt)
    # this default will not be used for TheiaCoV as that workflow series pass in the expected length based on the organism tag
    Int genome_length = 5000000 

    String? workflow_series

    # theiacov inputs
    Int? min_length
    Int? max_length
    String? run_prefix

    # kraken inputs
    String? target_organism
    Boolean call_kraken = false
    Int? kraken_disk_size
    Int? kraken_memory
    Int? kraken_cpu
    File? kraken_db

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
    if ((call_kraken) && defined(kraken_db)) {
      call kraken2.kraken2_standalone as kraken2_se {
        input:
          samplename = samplename,
          read1 = read1,
          kraken2_db = select_first([kraken_db]),
          disk_size = kraken_disk_size,
          memory = kraken_memory,
          cpu = kraken_cpu
      }
      call kraken2.kraken2_parse_classified as kraken2_recalculate_abundances {
        input:
          samplename = samplename,
          kraken2_report = kraken2_se.kraken2_report,
          kraken2_classified_report = kraken2_se.kraken2_classified_report
      }
    } if ((call_kraken) && ! defined(kraken_db)) {
      String kraken_db_warning = "Kraken database not defined"
    }

    # rasusa for random downsampling
    call rasusa_task.rasusa {
      input:
        read1 = read1,
        samplename = samplename,
        coverage = downsampling_coverage,
        genome_length = genome_length
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

  if (workflow_series == "theiaeuk") {
    # rasusa for random downsampling
    call rasusa_task.rasusa as theiaeuk_rasusa {
      input:
        read1 = read1,
        samplename = samplename,
        coverage = downsampling_coverage,
        genome_length = genome_length
    }

    # nanoq for filtering
    call nanoq_task.nanoq as theiaeuk_nanoq {
      input:
        read1 = theiaeuk_rasusa.read1_subsampled,
        samplename = samplename
    }
  }
  output { 
    # theiacov outputs
    # ncbi scrub outputs
    File? read1_dehosted = ncbi_scrub_se.read1_dehosted
    
    # kraken2 - theiacov and theiapro
    String kraken_version = select_first([kraken2_raw.version, kraken2_se.kraken2_version, ""])
    String kraken_docker = select_first([kraken2_raw.docker, kraken2_se.kraken2_docker, ""])
    Float? kraken_human = kraken2_recalculate_abundances_raw.percent_human
    Float? kraken_sc2 = kraken2_recalculate_abundances_raw.percent_sc2
    String? kraken_target_organism = kraken2_recalculate_abundances_raw.percent_target_organism
    String? kraken_target_organism_name = kraken2_raw.kraken_target_organism
    String kraken_report = select_first([kraken2_recalculate_abundances_raw.kraken_report, kraken2_recalculate_abundances.kraken_report, ""])
    Float? kraken_human_dehosted = kraken2_recalculate_abundances_dehosted.percent_human
    Float? kraken_sc2_dehosted = kraken2_recalculate_abundances_dehosted.percent_sc2
    String? kraken_target_organism_dehosted = kraken2_recalculate_abundances_dehosted.percent_target_organism
    File? kraken_report_dehosted = kraken2_recalculate_abundances_dehosted.kraken_report
    String kraken_database = select_first([kraken2_raw.database, kraken2_se.kraken2_database, kraken_db_warning, ""])
   
    # estimated genome length -- by default for TheiaProk this is 5Mb
    Int est_genome_length = genome_length

    # nanoq outputs
    File read1_clean = select_first([nanoq.filtered_read1, read_filtering.read1_clean, theiaeuk_nanoq.filtered_read1])
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