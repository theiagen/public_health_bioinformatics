version 1.0

task fetch_bs {
  input {
    String sample_name
    String basespace_sample_name
    String? basespace_sample_id
    String basespace_collection_id
    String api_server
    String access_token
    
    Int mem_size_gb = 8
    Int cpu = 2
    Int disk_size = 100
    Int preemptible = 1

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    # set basespace name and id variables
    if [[ ! -z "~{basespace_sample_id}" ]]; then
      sample_identifier="~{basespace_sample_name}"
      dataset_name="~{basespace_sample_id}"
    else
      sample_identifier="~{basespace_sample_name}"
      dataset_name="~{basespace_sample_name}"
    fi
    
    # print all relevant input variables to stdout
    echo -e "sample_identifier: ${sample_identifier}\ndataset_name: ${dataset_name}\nbasespace_collection_id: ~{basespace_collection_id}"
      
    #Set BaseSpace comand prefix
    bs_command="bs --api-server=~{api_server} --access-token=~{access_token}"
    echo "bs_command: ${bs_command}"

    #Grab BaseSpace Run_ID from given BaseSpace Run Name
    run_id=$(${bs_command} list run | grep "~{basespace_collection_id}" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
    echo "run_id: ${run_id}" 
    if [[ ! -z "${run_id}" ]]; then 
      #Grab BaseSpace Dataset ID from dataset lists within given run 
      dataset_id_array=($(${bs_command} list dataset --input-run=${run_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' )) 
      echo "dataset_id: ${dataset_id_array[*]}"
    else 
      #Try Grabbing BaseSpace Dataset ID from project name
      echo "Could not locate a run_id via Basespace runs, attempting to search Basespace projects now..."
      project_id=$(${bs_command} list project | grep "~{basespace_collection_id}" | awk -F "|" '{ print $3 }' | awk '{$1=$1;print}' )
      echo "project_id: ${project_id}" 
      if [[ ! -z "${project_id}" ]]; then 
        echo "project_id identified via Basespace, now searching for dataset_id within project_id ${project_id}..."
        dataset_id_array=($(${bs_command} list dataset --project-id=${run_id} | grep "${dataset_name}" | awk -F "|" '{ print $3 }' )) 
        echo "dataset_id: ${dataset_id_array[*]}"
      else       
        echo "No run or project id found associated with input basespace_collection_id: ~{basespace_collection_id}" >&2
        exit 1
      fi      
    fi

    #Download reads by dataset ID
    for index in ${!dataset_id_array[@]}; do
      dataset_id=${dataset_id_array[$index]}
      mkdir ./dataset_${dataset_id} && cd ./dataset_${dataset_id}
      echo "dataset download: ${bs_command} download dataset -i ${dataset_id} -o . --retry"
      ${bs_command} download dataset -i ${dataset_id} -o . --retry && cd ..
      echo -e "downloaded data: \n $(ls ./dataset_*/*)"
    done

    # rename FASTQ files to add back in underscores that Illumina/Basespace changed into hyphens
    echo "Concatenating and renaming FASTQ files to add back underscores in basespace_sample_name"
    # setting a new bash variable to use for renaming during concatenation of FASTQs
    SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES=$(echo $sample_identifier | sed 's|_|-|g' | sed 's|\.|-|g')

    #Combine non-empty read files into single file without BaseSpace filename cruft
    ##FWD Read
    lane_count=0
    for fwd_read in ./dataset_*/${SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES}_*R1_*.fastq.gz; do
      if [[ -s $fwd_read ]]; then
        echo "cat fwd reads: cat $fwd_read >> ~{sample_name}_R1.fastq.gz" 
        cat $fwd_read >> ~{sample_name}_R1.fastq.gz
        lane_count=$((lane_count+1))
      fi
    done
    ##REV Read
    for rev_read in ./dataset_*/${SAMPLENAME_HYPHEN_INSTEAD_OF_UNDERSCORES}_*R2_*.fastq.gz; do
      if [[ -s $rev_read ]]; then 
        echo "cat rev reads: cat $rev_read >> ~{sample_name}_R2.fastq.gz" 
        cat $rev_read >> ~{sample_name}_R2.fastq.gz
      fi
    done
    echo "Lane Count: ${lane_count}"
  >>>
  output {
    File read1 = "~{sample_name}_R1.fastq.gz"
    File? read2 = "~{sample_name}_R2.fastq.gz"
  }
  runtime {
    docker: docker
    memory: "~{mem_size_gb} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: preemptible
  }
}