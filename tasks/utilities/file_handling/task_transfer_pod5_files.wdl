version 1.0

task transfer_pod5_files {
  input {
    String pod5_bucket_path             # Terra bucket path (e.g., "gs://your-terra-bucket/pod5_uploads/")
    Int disk_size = 100
    Int memory = 32
    Int cpu = 8
    String docker =  "us-docker.pkg.dev/general-theiagen/cloudsdktool/google-cloud-cli:427.0.0-alpine"
  }
  command <<<
     set -euo pipefail

     # Create a directory for downloaded `.pod5` files
    mkdir -p pod5_downloads

    echo "Listing and downloading .pod5 files from ~{pod5_bucket_path}"
    gcloud storage ls -r "~{pod5_bucket_path}" | grep "\.pod5$" > pod5_files_list.txt

    # Check if any files are found
    if [ ! -s pod5_files_list.txt ]; then
      echo "ERROR: No POD5 files found in ~{pod5_bucket_path}" >&2
      exit 1
    fi

    # Download all `.pod5` files locally
    while read -r file_path; do
      local_path="pod5_downloads/$(basename "$file_path")"
      gcloud storage cp "$file_path" "$local_path" || { echo "ERROR: Failed to download $file_path"; exit 1; }
      echo "$local_path" >> downloaded_pod5_files.txt
    done < pod5_files_list.txt
  >>>

  output {
    Array[File] pod5_file_paths = read_lines("downloaded_pod5_files.txt")  # Local paths of downloaded `.pod5` files
  }
   runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 1
  }
}
