version 1.0

import "../../tasks/quality_control/read_filtering/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/contamination/task_kraken2.wdl" as kraken2
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow read_QC_trim_ont {
  meta {
    description: "Runs basic QC on Oxford Nanopore (ONT) reads with nanoplot, rasusa downsampling, and nanoq filtering"
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
    
    # ncbi-scrub inputs
    Int? ncbi_scrub_cpu
    Int? ncbi_scrub_disk_size
    String? ncbi_scrub_docker
    Int? ncbi_scrub_memory

    # artic guppyplex inputs
    Int? artic_guppyplex_cpu
    Int? artic_guppyplex_disk_size
    String? artic_guppyplex_docker
    Int? artic_guppyplex_memory

    # kraken inputs
    String? target_organism
    Boolean call_kraken = false
    Int? kraken_cpu
    File? kraken_db
    Int? kraken_disk_size
    String? kraken_docker_image
    Int? kraken_memory

    # parse classified reads inputs
    Int? kraken2_recalculate_abundances_cpu
    Int? kraken2_recalculate_abundances_disk_size
    String? kraken2_recalculate_abundances_docker
    Int? kraken2_recalculate_abundances_memory
  
    # rasusa downsampling inputs
    Float downsampling_coverage = 150
    Int? rasusa_cpu
    Int? rasusa_disk_size
    String? rasusa_docker
    Int? rasusa_memory
    String? rasusa_bases
    Int? rasusa_seed
    Float? rasusa_fraction_of_reads
    Int? rasusa_number_of_reads

    # nanoq inputs
    Int? nanoq_cpu
    Int? nanoq_disk_size
    String? nanoq_docker
    Int? nanoq_memory
    Int? nanoq_max_read_length
    Int? nanoq_min_read_length
    Int? nanoq_max_read_qual
    Int? nanoq_min_read_qual  
  }
  if ("~{workflow_series}" == "theiacov") {
    call ncbi_scrub.ncbi_scrub_se {
      input:
        read1 = read1,
        samplename = samplename,
        cpu = ncbi_scrub_cpu,
        disk_size = ncbi_scrub_disk_size,
        docker = ncbi_scrub_docker,
        memory = ncbi_scrub_memory
    }
    call artic_guppyplex.read_filtering {
      input:
        read1 = ncbi_scrub_se.read1_dehosted,
        samplename = samplename,
        min_length = min_length,
        max_length = max_length,
        run_prefix = run_prefix,
        cpu = artic_guppyplex_cpu,
        disk_size = artic_guppyplex_disk_size,
        docker = artic_guppyplex_docker,
        memory = artic_guppyplex_memory
    }
    call kraken2.kraken2_theiacov as kraken2_raw {
      input:
        samplename = samplename,
        read1 = read1,
        target_organism = target_organism,
        kraken2_db = kraken_db,
        disk_size = kraken_disk_size,
        memory = kraken_memory,
        cpu = kraken_cpu,
        docker_image = kraken_docker_image
    }
    call kraken2.kraken2_parse_classified as kraken2_recalculate_abundances_raw {
      input:
        samplename = samplename,
        kraken2_report = kraken2_raw.kraken_report,
        kraken2_classified_report = kraken2_raw.kraken2_classified_report,
        target_organism = target_organism,
        disk_size = kraken2_recalculate_abundances_disk_size,
        memory = kraken2_recalculate_abundances_memory,
        cpu = kraken2_recalculate_abundances_cpu,
        docker = kraken2_recalculate_abundances_docker
    }  
    call kraken2.kraken2_theiacov as kraken2_dehosted {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_se.read1_dehosted,
        target_organism = target_organism,
        kraken2_db = kraken_db,
        disk_size = kraken_disk_size,
        memory = kraken_memory,
        cpu = kraken_cpu,
        docker_image = kraken_docker_image
    }
    call kraken2.kraken2_parse_classified as kraken2_recalculate_abundances_dehosted {
      input:
        samplename = samplename,
        kraken2_report = kraken2_dehosted.kraken_report,
        kraken2_classified_report = kraken2_dehosted.kraken2_classified_report,
        target_organism = target_organism,
        disk_size = kraken2_recalculate_abundances_disk_size,
        memory = kraken2_recalculate_abundances_memory,
        cpu = kraken2_recalculate_abundances_cpu,
        docker = kraken2_recalculate_abundances_docker
    } 
  }
  if ("~{workflow_series}" == "theiaprok" || "~{workflow_series}" == "theiaeuk") {
    # rasusa for random downsampling
    call rasusa_task.rasusa {
      input:
        read1 = read1,
        samplename = samplename,
        bases = rasusa_bases,
        coverage = downsampling_coverage,
        cpu = rasusa_cpu,
        disk_size = rasusa_disk_size,
        docker = rasusa_docker,
        frac = rasusa_fraction_of_reads,
        genome_length = genome_length,
        memory = rasusa_memory,
        num = rasusa_number_of_reads,
        seed = rasusa_seed
    }
    # nanoq for filtering
    call nanoq_task.nanoq {
      input:
        read1 = rasusa.read1_subsampled,
        samplename = samplename,
        cpu = nanoq_cpu,
        disk_size = nanoq_disk_size,
        docker = nanoq_docker,
        max_read_length = nanoq_max_read_length,
        min_read_length = nanoq_min_read_length,
        max_read_qual = nanoq_max_read_qual,
        min_read_qual = nanoq_min_read_qual,
        memory = nanoq_memory 
    }
    if ("~{workflow_series}" == "theiaprok") {
      if (call_kraken && defined(kraken_db)) {
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
      } 
    if ((call_kraken) && ! defined(kraken_db)) {
        String kraken_db_warning = "Kraken database not defined"
      }
    }
  }
  output { 
    # theiacov outputs
    # ncbi scrub outputs
    File? read1_dehosted = ncbi_scrub_se.read1_dehosted
    
    # kraken2 - theiacov and theiaprok
    String kraken_version = select_first([kraken2_raw.version, kraken2_se.kraken2_version, ""])
    String kraken_docker = select_first([kraken2_raw.docker, kraken2_se.kraken2_docker, ""])
    Float? kraken_human = kraken2_recalculate_abundances_raw.percent_human
    String? kraken_sc2 = kraken2_recalculate_abundances_raw.percent_sc2
    String? kraken_target_organism = kraken2_recalculate_abundances_raw.percent_target_organism
    String? kraken_target_organism_name = kraken2_raw.kraken_target_organism
    String kraken_report = select_first([kraken2_recalculate_abundances_raw.kraken_report, kraken2_recalculate_abundances.kraken_report, ""])
    Float? kraken_human_dehosted = kraken2_recalculate_abundances_dehosted.percent_human
    String? kraken_sc2_dehosted = kraken2_recalculate_abundances_dehosted.percent_sc2
    String? kraken_target_organism_dehosted = kraken2_recalculate_abundances_dehosted.percent_target_organism
    File? kraken_report_dehosted = kraken2_recalculate_abundances_dehosted.kraken_report
    String kraken_database = select_first([kraken2_raw.database, kraken2_se.kraken2_database, kraken_db_warning, ""])
   
    # estimated genome length -- by default for TheiaProk this is 5Mb
    Int est_genome_length = genome_length

    # nanoq outputs
    File read1_clean = select_first([nanoq.filtered_read1, read_filtering.read1_clean])
    String? nanoq_version = nanoq.version

    # rasusa outputs
    String? rasusa_version = rasusa.rasusa_version  
  }
}