version 1.0

task list_files_by_extension {
  input {
    String bucket_path # GCS bucket path containing files (e.g., "gs://your-terra-bucket/uploads/")
    String file_extension # File extension to search for (e.g., ".pod5", ".fastq.gz")
    Int disk_size = 100
    Int memory = 4
    Int cpu = 1
    String docker =  "us-docker.pkg.dev/general-theiagen/cloudsdktool/google-cloud-cli:427.0.0-alpine"
  }
  command <<< 
    set -euo pipefail

    echo "Listing files with extension ~{file_extension} in ~{bucket_path}"
    gcloud storage ls -r "~{bucket_path}" | grep "~{file_extension}$" > files_list.txt

    # Check if any files are found
    if [ ! -s files_list.txt ]; then
      echo "ERROR: No files with extension ~{file_extension} found in ~{bucket_path}" >&2
      exit 1
    fi
  >>>
  output {
    Array[File] file_paths = read_lines("files_list.txt")
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
