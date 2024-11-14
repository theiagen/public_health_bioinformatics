version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task
import "../../tasks/basecalling/task_samtools_convert.wdl" as samtools_convert_task
import "../../tasks/basecalling/task_dorado_demux.wdl" as dorado_demux_task
import "../../tasks/utilities/file_handling/task_transfer_files.wdl" as transfer_fastq_files
import "../../tasks/utilities/data_import/task_create_terra_table.wdl" as terra_fastq_table
import "../../tasks/task_versioning.wdl" as versioning_task

workflow dorado_basecalling_workflow {
  meta {
    description: "GPU-accelerated workflow for basecalling Oxford Nanopore POD5 files, generating SAM outputs and supporting downstream demultiplexing and FASTQ output."
  }
  input {
    Array[File] input_files
    String dorado_model = "sup" 
    String kit_name
    String new_table_name
    String fastq_upload_path
    Boolean paired_end = false
    Boolean assembly_data = false
    String? file_ending
    String terra_project
    String terra_workspace
    String fastq_file_name             
  }
  call versioning_task.version_capture {
    input:
  }
  scatter (file in input_files) {
    call basecall_task.basecall as dorado_basecall {
      input:
        input_file = file,
        dorado_model = dorado_model,
        kit_name = kit_name,
    }
  }
  call samtools_convert_task.samtools_convert {
    input:
      sam_files = flatten(dorado_basecall.sam_files)
  }
  call dorado_demux_task.dorado_demux {
    input:
      bam_files = samtools_convert.bam_files,
      kit_name = kit_name,
      fastq_file_name = fastq_file_name,
      dorado_model_used = dorado_basecall.dorado_model_used
  }
  call transfer_fastq_files.transfer_files as transfer_files {
    input:
      files_to_transfer = dorado_demux.fastq_files,
      target_bucket = fastq_upload_path
  }
  if (defined(transfer_files.transferred_files)) {
    call terra_fastq_table.create_terra_table as create_terra_table {
      input:
        new_table_name = new_table_name,
        data_location_path = fastq_upload_path,
        paired_end = paired_end,
        assembly_data = assembly_data,
        file_ending = file_ending,
        terra_project = terra_project,
        terra_workspace = terra_workspace
    }
  }
  output {
    # Version Captures
    String dorado_phb_version = version_capture.phb_version
    String dorado_analysis_date = version_capture.date

    # Outputs from the main steps
    Array[File] fastq_files = dorado_demux.fastq_files
    File? terra_table_tsv = create_terra_table.terra_table_to_upload

    # Versions and model used 
    String dorado_version = dorado_demux.dorado_version
    String samtools_version = samtools_convert.samtools_version
    String dorado_model_used = dorado_demux.dorado_model_used
  }
}
