version 1.0

task subset_pod5s {
  input {
    Array[File] pod5_files
    String pod5_bucket_path

    Int cpu = 1
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/pod5:0.3.23"
    Int memory = 4
  }
  command <<<
    pod5 --version | tee VERSION

    # rename the pod5_bucket path to match the localization path
    # gs://bucket_path
    localized_path=$(echo "~{pod5_bucket_path}" | cut -d'/' -f3-)

    pod5 view /cromwell_root/${localized_path}/ --include "read_id, channel" --output summary.tsv
    pod5 subset /cromwell_root/${localized_path}/ --summary summary.tsv --columns channel --output split_by_channel
  >>>
  output {
    Array[File] pod5s_by_channel = glob("split_by_channel/*pod5")
    File summary = "summary.tsv"
    String pod5_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 1
  }
}