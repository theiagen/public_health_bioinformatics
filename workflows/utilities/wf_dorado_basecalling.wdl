version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as dorado_basecall_task
import "../../tasks/basecalling/task_dorado_demux.wdl" as dorado_demux_task
import "../../tasks/basecalling/task_dorado_trim.wdl" as dorado_trim_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_import/task_create_terra_table_from_array.wdl" as fastq_upload_table
import "../../tasks/utilities/file_handling/task_find_files.wdl" as find_files_task


workflow dorado_basecalling {
  meta {
    description: "GPU-accelerated workflow for basecalling Oxford Nanopore POD5 files, generating BAM outputs and supporting downstream demultiplexing and FASTQ output."
  }
  input {
    String pod5_bucket_path  # GCS bucket path containing POD5 files
    String dorado_model = "sup" 
    String kit_name
    String new_table_name
    String terra_project
    String terra_workspace
    String output_file_prefix 
    Boolean demux_notrim = false
    File? custom_primers
  }
  call versioning.version_capture {
    input:
  }
  call find_files_task.find_files {
    input:
      bucket_path = pod5_bucket_path,
      file_extension = ".pod5"
  }
  scatter (file_path in find_files.file_paths) {
    call dorado_basecall_task.dorado_basecall {
      input:
        pod5_file = file_path,
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
  if (defined(custom_primers)) {
    call dorado_trim_task.dorado_trim {
      input:
        fastq_files = dorado_demux.fastq_files,
        custom_primers = select_first([custom_primers])
    }
  }
  call fastq_upload_table.create_table_from_array {
    input:
      new_table_name = new_table_name,
      file_paths = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files]),
      file_ending = ".fastq.gz",
      output_file_column_name = "read1",
      data_source = "Dorado_Basecalling_PHB",
      terra_project = terra_project,
      terra_workspace = terra_workspace
  }
  output {
    String dorado_phb_version = version_capture.phb_version
    String dorado_analysis_date = version_capture.date
    Array[File] fastq_files = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files])
    File terra_table_tsv = create_table_from_array.terra_table_to_upload
    String dorado_version = dorado_demux.dorado_version
    String dorado_model_used = dorado_demux.dorado_model_name
  }
}
