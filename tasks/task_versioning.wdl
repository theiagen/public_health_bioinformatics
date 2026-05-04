version 1.0

task version_capture {
  input {
    String? timezone
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/alpine-plus-bash:3.20.0"

    String? workflow_name
  }
  meta {
    volatile: true
  }
  command {
    VERSION_TAG="v4.1.0"
    echo "PHB ${VERSION_TAG}" > PHB_VERSION

    export TZ=~{timezone}
    date -I > TODAY
    
    if [ -z "~{workflow_name}" ]; then
      echo "No workflow name provided, skipping default input capture."

    else
      echo "getting default inputs for the current workflow from the PHB documentation"
      wget https://raw.githubusercontent.com/theiagen/public_health_bioinformatics/refs/tags/${VERSION_TAG}/docs/assets/tables/all_inputs.tsv
    
      cut -f-2,5,7-8 all_inputs.tsv | grep ~{workflow_name} > workflow_inputs.tsv 
      awk -F'\t' '$3 != ""' workflow_inputs.tsv > default_workflow_inputs.tsv

      echo "# DEFAULT DATABASES FOR ~{workflow_name}" > default_inputs.tsv
      grep "database" default_workflow_inputs.tsv | cut -f-3 >> default_inputs.tsv
      
      echo -e "\n# DEFAULT REFERENCE FILES FOR ~{workflow_name}" >> default_inputs.tsv
      grep "reference" default_workflow_inputs.tsv | cut -f-3 >> default_inputs.tsv
      
      echo -e "\n# DEFAULT DOCKER IMAGES FOR ~{workflow_name}" >> default_inputs.tsv
      grep "docker" default_workflow_inputs.tsv | cut -f-3 >> default_inputs.tsv
      
    fi
  }
  output {
    String date = read_string("TODAY")
    String phb_version = read_string("PHB_VERSION")
    
    File? default_workflow_inputs = "default_inputs.tsv"
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

