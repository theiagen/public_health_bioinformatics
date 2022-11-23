version 1.0

task version_capture {
  input {
    String? timezone
  }
  meta {
    volatile: true
  }
  command {
    PHBG_Version="PHBG v1.0.0"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$PHBG_Version" > PHBG_VERSION
  }
  output {
    String date = read_string("TODAY")
    String phbg_version = read_string("PHBG_VERSION")
  }
  runtime {
    memory: "1 GB"
    cpu: 1
    docker: "quay.io/theiagen/utility:1.1"
    disks: "local-disk 10 HDD"
    dx_instance_type: "mem1_ssd1_v2_x2" 
  }
}
