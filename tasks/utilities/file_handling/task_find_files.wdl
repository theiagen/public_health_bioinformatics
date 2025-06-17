version 1.0

task find_files {
  input {
    String bucket_path # include trailing slash if directory
    String file_extension

    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/cloudsdktool/google-cloud-cli:427.0.0-alpine"
    Int memory = 2
  }
  command <<< 
    set -euo pipefail

    echo "Listing files with extension ~{file_extension} in ~{bucket_path}"
    gcloud storage ls -r "~{bucket_path}" | grep "~{file_extension}$" > file_list.txt

    if [ ! -s file_list.txt ]; then
      echo "ERROR: No files with extension ~{file_extension} found in ~{bucket_path}" >&2
      exit 1
    fi
  >>>
  output {
    Array[String] file_paths = read_lines("file_list.txt")
  }
   runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 3
  }
}
