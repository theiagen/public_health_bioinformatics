version 1.0

task create_table_from_array {
  input {
    String new_table_name
    Array[File] file_paths
    String file_ending
    String output_file_column_name
    String? data_source

    Map[String, String?] columns_to_export

    String terra_project
    String terra_workspace

    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21"
    Int memory = 2
  }
    File columns_to_export_json = write_json(columns_to_export)
  meta {
    volatile: true
  }
  String new_table_name_updated = sub(new_table_name, " ", "_")
  command <<<
    echo "DEBUG: starting to creating a terra table"

    # writing these to stderr so if errors occur, the user can check what they provided.
    echo "DEBUG: new_table_name: ~{new_table_name_updated}" >&2
    echo "DEBUG: file_ending: ~{file_ending}" >&2
    echo "DEBUG: file_column_name: ~{output_file_column_name}" >&2
    echo "DEBUG: data_source: ~{data_source}" >&2
    
    # add additional columns to the terra table
    jq -r '[.[] | .left], [.[] | .right] | @tsv' ~{columns_to_export_json} > exported_columns.tsv

    column_names=$(head -n1 exported_columns.tsv)
    column_content=$(tail -n1 exported_columns.tsv)
 
    echo -e "entity:~{new_table_name_updated}_id\t~{output_file_column_name}\tupload_date\ttable_created_by\t${column_names}" > terra_table_to_upload.tsv

    UPLOAD_DATE=$(date -I)

    # loop through array of files and extract their basename to add to the terra table
    for file in ~{sep=" " file_paths}; do
      file_basename=$(basename "$file" ~{file_ending})

      echo -e "${file_basename}\t${file//*\/cromwell_root/gs:\/}\t${UPLOAD_DATE}\t~{data_source}\t${column_content}" >> terra_table_to_upload.tsv
    done
    
    echo "DEBUG: terra_table_to_upload.tsv created"
    # set error handling to exit if the subsequent import_large_tsv.py task fails
    set -euo pipefail

    python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --tsv terra_table_to_upload.tsv
  >>>
  output {
    File terra_table_to_upload = "terra_table_to_upload.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}