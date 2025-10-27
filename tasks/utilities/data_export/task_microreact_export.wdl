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

    python /scripts/microreact_export.py \
      --project_name ~{project_name} \
      ~{"--project_url " + project_url} \
      --metadata_tsv ~{metadata_tsv} \
      --id_column ~{id_column} \
      ~{"--date_column " + date_column} \
      --tree_files ~{sep=" " tree_files} \
      --selected_columns ~{sep=" " metadata_columns} \
      ~{"--access_token " + access_token} \
      --restricted_access ~{restricted_access} \
      --update ~{update_project} \
      --remove_file_columns ~{remove_file_columns}
      
  >>>
  output {
    File microreact_json = "project_input.json"
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