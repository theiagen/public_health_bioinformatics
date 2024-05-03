version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
  }
  meta {
    volatile: true
  }
  command {
    PHB_Version="PHB v2.0.1: branch smw-tb-2024-05-03-dev"
    ~{default='' 'export TZ=' + timezone}
    date +"%Y-%m-%d" > TODAY
    echo "$PHB_Version" > PHB_VERSION
  }
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
  }
}

