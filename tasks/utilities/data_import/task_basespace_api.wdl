version 1.0

task fetch_bs {
  input {
    String sample_name
    String basespace_sample_name
    String basespace_collection_id
    String access_token

    Boolean validate_paired_end = true
    Boolean validate_lane_naming = false
    Boolean group_by_lane = false

    Int memory = 8
    Int cpu = 2
    Int disk_size = 100

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bioforklift:0.5.0-dev"
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    from bioforklift.basespace import BaseSpace

    # convert boolean WDL inputs to pythonic type bool
    validate_paired_end = "~{validate_paired_end}" == "true"
    validate_lane_naming = "~{validate_lane_naming}" == "true"
    group_by_lane = "~{group_by_lane}" == "true"

    bs = BaseSpace(
        access_token="~{access_token}",
    )

    bs.fetch_sample_fastqs(
        collection_id="~{basespace_collection_id}",
        samples=["~{basespace_sample_name}"],
        dataset_types=["common.fastq"],             # required for this task
        concatenate=True,                           # required for this task
        group_by_lane=group_by_lane,
        validate_paired_end=validate_paired_end,
        validate_lane_naming=validate_lane_naming,
    )
    CODE

    # Cannot rename the concatenated FASTQs if they have the same name
    if [[ "~{basespace_sample_name}" != "~{sample_name}" ]]; then
        mv "~{basespace_sample_name}_R1.fastq.gz" "~{sample_name}_R1.fastq.gz"
        mv "~{basespace_sample_name}_R2.fastq.gz" "~{sample_name}_R2.fastq.gz"
    fi
  >>>
  output {
    File read1 = "~{sample_name}_R1.fastq.gz"
    File read2 = "~{sample_name}_R2.fastq.gz"
    File basespace_log = "bioforklift.log"
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
