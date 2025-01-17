version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as basecall_task
import "../../tasks/basecalling/task_dorado_demux.wdl" as dorado_demux_task
import "../../tasks/utilities/file_handling/task_transfer_files.wdl" as transfer_fastq_files
import "../../tasks/utilities/data_import/task_create_terra_table.wdl" as terra_fastq_table
import "../../tasks/task_versioning.wdl" as versioning_task
import "../../tasks/utilities/file_handling/task_list_pod5_files.wdl" as list_pod5_files_task    
import "../../tasks/basecalling/task_dorado_trim.wdl" as task_dorado_trim     

workflow dorado_basecalling_workflow {
  meta {
    description: "GPU-accelerated workflow for basecalling Oxford Nanopore POD5 files, generating BAM outputs and supporting downstream demultiplexing and FASTQ output."
  }
  input {
    String pod5_bucket_path  # GCS bucket path containing POD5 files
    String dorado_model = "sup" 
    String kit_name
    String new_table_name
    String fastq_upload_path
    String terra_project
    String terra_workspace
    String output_file_prefix 
    Boolean demux_notrim = false
    File? custom_primers
  }

  call versioning_task.version_capture {
    input:
  }

  call list_pod5_files_task.list_pod5_files as list_pod5 {
    input:
      pod5_bucket_path = pod5_bucket_path
  }

  scatter (pod5_path in list_pod5.pod5_file_paths) {
    call basecall_task.basecall as dorado_basecall {
      input:
        input_file = pod5_path,
        dorado_model = dorado_model,
        kit_name = kit_name
    }
  }

  call dorado_demux_task.dorado_demux {
    input:
      bam_files = flatten(dorado_basecall.bam_files),
      kit_name = kit_name,
      output_file_prefix = output_file_prefix,
      dorado_model_used = dorado_basecall.dorado_model_used[0],
      demux_notrim = demux_notrim
  }

  # Optional trimming step
  if (defined(custom_primers)) {
    call task_dorado_trim.dorado_trim as dorado_trim {
      input:
        fastq_files = dorado_demux.fastq_files,
        custom_primers = select_first([custom_primers])
    }
  }
  call transfer_fastq_files.transfer_files {
    input:
      files_to_transfer = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files]),
      target_bucket = fastq_upload_path
  }

  if (defined(transfer_files.transferred_files)) {
    call terra_fastq_table.create_terra_table {
      input:
        new_table_name = new_table_name,
        data_location_path = fastq_upload_path,
        paired_end = false,                  
        assembly_data = false,  
        file_ending = ".fastq.gz",         
        terra_project = terra_project,
        terra_workspace = terra_workspace
    }
  }
  output {
    # Version Captures
    String dorado_phb_version = version_capture.phb_version
    String dorado_analysis_date = version_capture.date

    # Outputs from the main steps
    Array[File] fastq_files = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files])
    File? terra_table_tsv = create_terra_table.terra_table_to_upload

    # Versions and model used 
    String dorado_version = dorado_demux.dorado_version
    String dorado_model_used = dorado_demux.dorado_model_name
  }
}
