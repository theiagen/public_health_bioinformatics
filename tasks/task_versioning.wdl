version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"
    String logging_api_url = "https://wdl-logger-ki7cuiofna-uc.a.run.app/log-task"  # Replace with actual URL
    String task_id = "version_capture_1.0"  # Unique Task ID
  }
  meta {
    volatile: true
  }
  command <<<
    PHB_Version="PHB v2.3.0"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$PHB_Version" > PHB_VERSION

    # Send log to Cloud Run API
    curl -X POST "${logging_api_url}" \
      -H "Content-Type: application/json" \
      -d '{
        "task_id": "'${task_id}'",
        "timestamp": "'$(date -u +%Y-%m-%dT%H:%M:%SZ)'",
        "status": "completed",
        "metadata": {
          "date": "'$(cat TODAY)'",
          "phb_version": "'$(cat PHB_VERSION)'"
        }
      }'
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
