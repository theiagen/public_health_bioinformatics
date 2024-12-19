version 1.0

task list_pod5_files {
  input {
    String pod5_bucket_path             # GCS bucket path containing `.pod5` files (e.g., "gs://your-terra-bucket/pod5_uploads/")
    Int disk_size = 100
    Int memory = 4
    Int cpu = 1
    String docker =  "us-docker.pkg.dev/general-theiagen/cloudsdktool/google-cloud-cli:427.0.0-alpine"
  }
  command <<<
    set -euo pipefail

    echo "Listing .pod5 files in ~{pod5_bucket_path}"
    gcloud storage ls -r "~{pod5_bucket_path}" | grep "\.pod5$" > pod5_files_list.txt

    # Check if any files are found
    if [ ! -s pod5_files_list.txt ]; then
      echo "ERROR: No POD5 files found in ~{pod5_bucket_path}" >&2
      exit 1
    fi
  >>>
  output {
    Array[File] pod5_file_paths = read_lines("pod5_files_list.txt")
  }
   runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 1
    maxRetries: 1
  }
}
