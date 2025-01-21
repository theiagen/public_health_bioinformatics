version 1.0

task check_gcp_project {
    input {
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/cloud-sdk:dev"
    }
    command <<<
        bash /scripts/check_project.sh > gcp_project.txt
    >>>
    output {
        File gcp_project = "gcp_project.txt"
    }
    runtime {
        docker: docker
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}