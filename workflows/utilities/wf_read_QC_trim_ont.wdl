version 1.0

import "../../tasks/quality_control/read_filtering/task_artic_guppyplex.wdl" as artic_guppyplex
import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as ete4_taxon_id
import "../../tasks/quality_control/read_filtering/task_nanoq.wdl" as nanoq_task
import "../../tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as ncbi_scrub
import "../../tasks/taxon_id/contamination/task_metabuli.wdl" as metabuli
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
    Int? metabuli_cpu
    File? metabuli_db
    Int? metabuli_disk_size
    String? metabuli_docker_image
    Int? metabuli_memory
    File? metabuli_taxdump_path = "gs://theiagen-public-resources-rp/reference_data/databases/metabuli/ncbi_taxdump_20260211.tar.gz"

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
  if (defined(target_organism)) {
    if (select_first([target_organism]) != "") {
      call ete4_taxon_id.ete4_taxon_id {
        input:
          taxon = select_first([target_organism]),
      }
    }
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
    call metabuli.metabuli as metabuli_theiacov_raw {
      input:
        samplename = samplename,
        read1 = read1,
        metabuli_db = metabuli_db,
        seq_mode = 3,
        taxon_id = ete4_taxon_id.taxon_id,
        disk_size = metabuli_disk_size,
        memory = metabuli_memory,
        cpu = metabuli_cpu,
        docker = metabuli_docker_image,
    }
    call metabuli.metabuli as metabuli_theiacov_dehosted {
      input:
        samplename = samplename,
        read1 = ncbi_scrub_se.read1_dehosted,
        taxon_id = ete4_taxon_id.taxon_id,
        metabuli_db = metabuli_db,
        seq_mode = 3,
        disk_size = metabuli_disk_size,
        memory = metabuli_memory,
        cpu = metabuli_cpu,
        docker = metabuli_docker_image,
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
    if ("~{workflow_series}" == "theiaprok" && defined(metabuli_db)) {
      call metabuli.metabuli as metabuli_theiaprok {
        input:
          samplename = samplename,
          read1 = read1,
          metabuli_db = metabuli_db,
          seq_mode = 3,
          taxon_id = ete4_taxon_id.taxon_id,
          disk_size = metabuli_disk_size,
          memory = metabuli_memory,
          cpu = metabuli_cpu,
      } 
    }
  }
  output { 
    # theiacov outputs
    # ncbi scrub outputs
    File? read1_dehosted = ncbi_scrub_se.read1_dehosted
    
    # kraken2 - theiacov and theiaprok
    String metabuli_version = select_first([metabuli_theiacov_raw.metabuli_version, metabuli_theiaprok.metabuli_version, ""])
    String metabuli_docker = select_first([metabuli_theiacov_raw.metabuli_docker, metabuli_theiaprok.metabuli_docker, ""])
    Float? metabuli_percent_human = select_first([metabuli_theiacov_raw.metabuli_percent_human, metabuli_theiaprok.metabuli_percent_human, ""])
    String? metabuli_percent_sc2 = select_first([metabuli_theiacov_raw.metabuli_percent_sc2, metabuli_theiaprok.metabuli_percent_sc2, ""])
    String? metabuli_target_organism = ete4_taxon_id.taxon_name
    String? metabuli_percent_target_organism = select_first([metabuli_theiacov_raw.metabuli_percent_target_lineage, metabuli_theiaprok.metabuli_percent_target_lineage, ""])
    String? metabuli_taxon_id = ete4_taxon_id.taxon_id
    String metabuli_report = select_first([metabuli_theiacov_raw.metabuli_report, metabuli_theiaprok.metabuli_report, ""])
    Float? metabuli_percent_human_dehosted = metabuli_theiacov_dehosted.metabuli_percent_human
    String? metabuli_percent_sc2_dehosted = metabuli_theiacov_dehosted.metabuli_percent_sc2
    String? metabuli_percent_target_organism_dehosted = metabuli_theiacov_dehosted.metabuli_percent_target_lineage
    File? metabuli_report_dehosted = metabuli_theiacov_dehosted.metabuli_report
    String metabuli_database = select_first([metabuli_theiacov_raw.metabuli_database, metabuli_theiaprok.metabuli_database, ""])
   
    # estimated genome length -- by default for TheiaProk this is 5Mb
    Int est_genome_length = genome_length

    # nanoq outputs
    File read1_clean = select_first([nanoq.filtered_read1, read_filtering.read1_clean])
    String? nanoq_version = nanoq.version

    # rasusa outputs
    String? rasusa_version = rasusa.rasusa_version  
  }
}