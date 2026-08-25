version 1.0

task fetch_bs {
  input {
    String sample_name
    String basespace_sample_name
    String basespace_collection_id
    String basespace_access_token
    String basespace_api_url = "https://api.basespace.illumina.com"

    Boolean validate_paired_end = true
    Boolean group_by_lane = true

    Int memory = 8
    Int cpu = 2
    Int disk_size = 250

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bioforklift:0.5.2-dev"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    set -euo pipefail

    python3 <<CODE
    from bioforklift.basespace import BaseSpace

    # convert boolean WDL inputs to pythonic type bool
    validate_paired_end = "~{validate_paired_end}" == "true"
    group_by_lane = "~{group_by_lane}" == "true"

    bs = BaseSpace(
        access_token="~{basespace_access_token}",
        basespace_api_url="~{basespace_api_url}"
    )

    bs.fetch_sample_fastqs(
        collection_id="~{basespace_collection_id}",
        samples=["~{basespace_sample_name}"],
        validate_paired_end=validate_paired_end,
        group_by_lane=group_by_lane,
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
