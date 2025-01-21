version 1.0

task check_gcp_project {
    input {
        String docker = "us-docker.pkg.dev/general-theiagen/theiagen/cloud-sdk:dev"
    }
    command <<<
        project_id=$(gcloud config get-value project)
        echo "$project_id" > PROJECT_ID
    >>>
    output {
        String gcp_project = read_string("PROJECT_ID")
    }
    runtime {
        docker: docker
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}