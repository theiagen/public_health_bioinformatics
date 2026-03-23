version 1.0

import "../../tasks/taxon_id/task_ete4_taxon_id.wdl" as identify_taxon_id_task
import "../../tasks/utilities/task_datasets_genome_length.wdl" as genome_length_task
import "../../tasks/quality_control/basic_statistics/task_cg_pipeline.wdl" as cg_pipeline_task
import "../../tasks/utilities/task_rasusa.wdl" as rasusa_task

workflow lod_table_process {
  input {
    String taxon
    File read1
    File read2
    String samplename
    Float? downsampling_level
  }
  # get the taxon id
  call identify_taxon_id_task.ete4_taxon_id as ete4_identify {
    input:
      taxon = taxon,
  }
  # get average genome length for the taxon
  call genome_length_task.datasets_genome_length as est_genome_length {
    input:
      taxon = ete4_identify.raw_taxon_id,
      use_ncbi_virus = false,
      complete = true,
      refseq = true
  }
  # get estimated coverage for downsampling
  call cg_pipeline_task.cg_pipeline as cg_pipeline {
    input:
      read1 = read1,
      read2 = read2,
      samplename = samplename,
      genome_length = est_genome_length.avg_genome_length
  }
  # only downsample and process if at appropriate coverage level (ignore original sample)
  Boolean should_downsample = defined(downsampling_level) && select_first([downsampling_level, 0.0]) <= cg_pipeline.est_coverage
  if (should_downsample) {
    call rasusa_task.rasusa as rasusa {
      input:
        read1 = read1,
        read2 = read2,
        samplename = samplename,
        genome_length = est_genome_length.avg_genome_length,
        coverage = downsampling_level,
    }
    String rasusa_completed = "COMPLETED"
  }
  if (!should_downsample) {
    String rasusa_skipped = "SKIPPED"
  }
  output {
    # ete4 - taxon identification
    String ncbi_taxon_id = ete4_identify.taxon_id
    String ncbi_taxon_name = ete4_identify.taxon_name
    String ncbi_read_extraction_rank = ete4_identify.taxon_rank
    String ete4_status = ete4_identify.ete4_status
    # String? ete4_version = ete4_identify.ete4_version
    # String? ete4_docker = ete4_identify.ete4_docker
    # NCBI datasets genome length estimation
    String taxon_avg_genome_length = est_genome_length.avg_genome_length
    # String? datasets_genome_length_docker = est_genome_length.ncbi_datasets_docker
    # String? datasets_genome_length_version = est_genome_length.ncbi_datasets_version
    # Read QC - cg pipeline outputs
    Float cg_pipeline_est_coverage = cg_pipeline.est_coverage
    File cg_pipeline_report = cg_pipeline.cg_pipeline_report
    # String cg_pipeline_docker = cg_pipeline.cg_pipeline_docker
    # Float? r1_mean_q_raw = cg_pipeline_raw.r1_mean_q
    # Float? r2_mean_q_raw = cg_pipeline_raw.r2_mean_q
    # Float? r1_mean_readlength_raw = cg_pipeline_raw.r1_mean_readlength
    # Float? r2_mean_readlength_raw = cg_pipeline_raw.r2_mean_readlength
    # Rasusa - downsampling of reads
    String rasusa_status = select_first([rasusa_completed, rasusa_skipped])
    File? read1_subsampled = rasusa.read1_subsampled
    File? read2_subsampled = rasusa.read2_subsampled
    # String? rasusa_version = rasusa.rasusa_version
  }
}