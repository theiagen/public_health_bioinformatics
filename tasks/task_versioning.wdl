version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"
  }
  meta {
    volatile: true
  }
  command <<<
    # VERSION_TAG is managed manually only on version updates
    VERSION_TAG="v4.3.0"
    # BRANCH_TAG is managed by CI; do NOT edit manually
    BRANCH_TAG="tj-basespace-wdl-dev"
    if [ -n "${BRANCH_TAG}" ]; then
      echo "PHB ${VERSION_TAG}; branch: ${BRANCH_TAG}" > PHB_VERSION
    else
      echo "PHB ${VERSION_TAG}" > PHB_VERSION
    fi

    export TZ=~{timezone}
    date -I > TODAY
  >>>
  output {
    String date = read_string("TODAY")
    String phb_version = read_string("PHB_VERSION")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: docker
    disks: "local-disk 10 HDD"
    disk: "10 GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 1
  }
}
