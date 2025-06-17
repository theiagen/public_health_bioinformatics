version 1.0

import "../../tasks/basecalling/task_dorado_basecall.wdl" as dorado_basecall_task
import "../../tasks/basecalling/task_dorado_demux.wdl" as dorado_demux_task
import "../../tasks/basecalling/task_dorado_trim.wdl" as dorado_trim_task
import "../../tasks/task_versioning.wdl" as versioning
import "../../tasks/utilities/data_import/task_array_to_terra.wdl" as create_fastq_table
import "../../tasks/utilities/file_handling/task_find_files.wdl" as find_files_task
import "../../tasks/utilities/file_handling/task_chunk_files.wdl" as chunk_file

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
    Boolean demux_no_trim = false
    File? custom_primers

    Int number_chunks = 4
  }
  call versioning.version_capture {
    input:
  }
  call find_files_task.find_files {
    input:
      bucket_path = pod5_bucket_path,
      file_extension = ".pod5"
  }
  Int number_of_files = length(find_files.file_paths)
  Int chunk_size = number_of_files / number_chunks
  if (number_chunks > number_of_files) {
    # if this condition is true, the scatter would be sad
    Int alt_number_chunks = number_of_files
  }
  scatter (i in range(select_first([alt_number_chunks, number_chunks]))) {
    Int start_line = (i * chunk_size + 1)
    Int end_line = if (i == number_chunks - 1) then number_of_files else (i + 1) * chunk_size
    # since wdl v1 doesn't do array splicing (:sob:), we have to mimic that behavior with line splicing out of a file 
    call chunk_file.chunk_files {
      input:
        file_list = write_lines(find_files.file_paths),
        start_line = start_line,
        end_line = end_line
    }
    call dorado_basecall_task.dorado_basecall {
      input:
        pod5_files = read_lines(chunk_files.chunk),
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
      demux_no_trim = demux_no_trim
  }
  if (defined(custom_primers)) {
    call dorado_trim_task.dorado_trim {
      input:
        fastq_files = dorado_demux.fastq_files,
        custom_primers = select_first([custom_primers])
    }
  }
  call create_fastq_table.create_table_from_array {
    input:
      new_table_name = new_table_name,
      file_paths = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files]),
      file_ending = ".fastq.gz",
      output_file_column_name = "read1",
      data_source = "Dorado_Basecalling_PHB",
      columns_to_export = {
        "dorado_docker": dorado_basecall.dorado_docker[0], 
        "dorado_version": dorado_basecall.dorado_version[0],
        "dorado_model_name": dorado_demux.dorado_model_name, 
        "dorado_basecalling_phb_verrsion": version_capture.phb_version, 
        "dorado_basecalling_analysis_date": version_capture.date
      },
      terra_project = terra_project,
      terra_workspace = terra_workspace
  }
  output {
    # final outputs
    Array[File] fastq_files = select_first([dorado_trim.trimmed_fastq_files, dorado_demux.fastq_files])
    String dorado_model_used = dorado_demux.dorado_model_name
    # task versioning
    String dorado_basecall_version = dorado_basecall.dorado_version[0]
    String dorado_basecall_docker = dorado_basecall.dorado_docker[0]
    String dorado_demux_version = dorado_demux.dorado_version
    String? dorado_trim_version = dorado_trim.dorado_version
    # uploaded table
    File terra_table_tsv = create_table_from_array.terra_table_to_upload
    # workflow versioning
    String dorado_basecalling_phb_version = version_capture.phb_version
    String dorado_basecalling_analysis_date = version_capture.date
  }
}
