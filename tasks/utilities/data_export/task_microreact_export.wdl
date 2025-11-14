version 1.0

task microreact_export {
  meta {
    description: "This task generates a Microreact input file from supplied metadata and tree files, submitting to Microreact when access token is provided."
  }
  input {
    String project_name
    String id_column
    File metadata_tsv
    String? date_column
    Array[String]? metadata_columns
    Array[File]? tree_files
    Array[File]? matrix_files
    String? access_token
    Boolean restricted_access = true
    Boolean update_project = false
    Boolean remove_file_columns = true
    String? project_url
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/microreact_export:0.1.0-dev"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 2
  }
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail
    
    tree_array=(~{sep=' ' tree_files})
    matrix_array=(~{sep=' ' matrix_files})
    metadata_column_array=(~{sep=' ' metadata_columns})

    python /scripts/microreact_export.py \
      --project_name ~{project_name} \
      ~{"--project_url " + project_url} \
      --metadata_tsv ~{metadata_tsv} \
      ~{if defined(matrix_files) && length(select_first([matrix_files, []])) > 0
        then "--matrix_files " else ""} "${matrix_array[@]}" \
      --id_column ~{id_column} \
      ~{"--date_column " + date_column} \
      ~{if defined(tree_files) && length(select_first([tree_files, []])) > 0 
        then "--tree_files " else ""} "${tree_array[@]}" \
      ~{if defined(metadata_columns) && length(select_first([metadata_columns, []])) > 0
          then "--selected_columns " else ""} "${metadata_column_array[@]}" \
      ~{if defined(access_token) then "--access_token " + access_token else ""} \
      ~{true="--restricted_access" false="" restricted_access} \
      ~{true="--update" false="" update_project} \
      ~{true="--remove_file_columns" false="" remove_file_columns}
      
  >>>
  output {
    File microreact_input = "~{project_name}_input.microreact"
    File? microreact_api_response = "microreact_response.json"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}