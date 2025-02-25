version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash-curl:3.20.0"
    String logging_api_url = "https://wdl-logger-932099266299.us-central1.run.app/log-task"
    String task_id = "version_capture_1.0"
    String workflow_name = "unknown"
  }
  meta {
    volatile: true
  }
command <<< 
    PHB_Version="PHB v2.3.0"

    # Ensure timezone is properly set (avoid empty variable issue)
    if [[ -n "~{timezone}" ]]; then
      export TZ="~{timezone}"
    fi

    date +"%Y-%m-%d" > TODAY
    echo "$PHB_Version" > PHB_VERSION

    IDENTITY_TOKEN=$(curl -s -X GET -H "Metadata-Flavor: Google" \
  "http://metadata/computeMetadata/v1/instance/service-accounts/default/identity?audience=${logging_api_url}")

    curl -X POST "~{logging_api_url}" \
          -H "Authorization: Bearer $IDENTITY_TOKEN" \
          -H "Content-Type: application/json" \
          -d '{ "task_id": "~{task_id}", "timestamp": "'$(date -u +%Y-%m-%dT%H:%M:%SZ)'", "workflow_name": "~{workflow_name}", "phb_version": "'$PHB_Version'" }'

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
    dx_instance_type: "mem1_ssd1_v2_x2" 
    preemptible: 1
  }
}
