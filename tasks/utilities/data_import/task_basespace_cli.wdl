version 1.0

task fetch_bs {
  input {
    String basespace_sample_name
    String? basespace_run_id
    String? basespace_project_id
    String access_token
    String api_server = "https://api.basespace.illumina.com"

    Int memory = 8
    Int cpu = 2
    Int disk_size = 100

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/basespace_cli:1.2.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    set -euo pipefail

    # set BaseSpace command prefix; --retry will retry transient errors up to 5 times
    bs_command="bs --api-server=~{api_server} --access-token=~{access_token} --retry"

    # exactly one of basespace_run_id / basespace_project_id must be provided
    if [[ -z "~{basespace_run_id}" && -z "~{basespace_project_id}" ]]; then
      echo "ERROR: must supply 'basespace_run_id' or 'basespace_project_id'" >&2
      exit 1
    fi
    if [[ -n "~{basespace_run_id}" && -n "~{basespace_project_id}" ]]; then
      echo "ERROR: supply only one of 'basespace_run_id' or 'basespace_project_id'" >&2
      exit 1
    fi

    # determine whether we're looking for datasets within a BaseSpace "Run" or "Project", and set appropriate variables for later use
    # runs use {{.ExperimentName}} (the user-given label); {{.Name}} is the instrument run id
    if [[ -n "~{basespace_run_id}" ]]; then
      bs_app="run"
      bs_app_id="~{basespace_run_id}"
      bs_dataset_flag="--input-run"
      bs_template_value="{{.ExperimentName}},{{.Id}}"
    else
      bs_app="project"
      bs_app_id="~{basespace_project_id}"
      bs_dataset_flag="--project-id"
      bs_template_value="{{.Name}},{{.Id}}"
    fi
    echo "Using BaseSpace ${bs_app} ID: '${bs_app_id}'"


    # attempt to resolve the run_id/project_id from the user input, which could be either the run_id, run name, project_id, or project name
    if result=$(${bs_command} get "${bs_app}" --name="${bs_app_id}" --template="${bs_template_value}" 2>/dev/null); then
      : # found by name
    elif result=$(${bs_command} get "${bs_app}" --id="${bs_app_id}" --template="${bs_template_value}" 2>/dev/null); then
      : # found by ID
    else
      echo "ERROR: Could not find '${bs_app_id}' in BaseSpace ${bs_app}" >&2
      exit 1
    fi

    #Download reads by dataset ID
    for index in ${!dataset_id_array[@]}; do
      dataset_id=${dataset_id_array[$index]}
      mkdir ./dataset_${dataset_id} && cd ./dataset_${dataset_id}
      echo "dataset download: ${bs_command} download dataset -i ${dataset_id} -o . --retry"
      ${bs_command} download dataset --retry -i ${dataset_id} -o . --retry && cd ..
      echo -e "downloaded data: \n $(ls ./dataset_*/*)"
    done
    # extract resolved name and ID from the result
    resolved_name=$(echo "${result}" | awk -F',' '{print $1}')
    resolved_id=$(echo "${result}" | awk -F',' '{print $2}')
    echo "Resolved BaseSpace ${bs_app} [name, id]: ['${resolved_name}', '${resolved_id}']"

    # rename FASTQ files to add back in underscores that Illumina/Basespace changed into hyphens
    echo "Concatenating and renaming FASTQ files to add back underscores in basespace_sample_name"
    # setting a new bash variable to use for renaming during concatenation of FASTQs
    for elm in ./dataset_${dataset_id}/*.fastq.gz; do
      echo "Checking Basespace file: $elm"
      filename=$(basename "$elm")

      if [[ "$filename" =~ [-] && "$sample_identifier" =~ [_] ]]; then
        echo "Basespace sample name for $filename contains dashes, input sample identifier $sample_identifier contains underscores, renaming identifier..."
        SAMPLENAME_RENAMED=$(echo "$sample_identifier" | sed 's|_|-|g' | sed 's|\.|-|g')
      fi
      if [[ "$filename" =~ [_] && "$sample_identifier" =~ [-] ]]; then
        echo "Basespace sample name for $filename contains underscores, input sample identifier $sample_identifier contains dashes, renaming identifier..."
        SAMPLENAME_RENAMED=$(echo "$sample_identifier" | sed 's|-|_|g')
      fi
      if [[ ("$filename" =~ [_] && "$sample_identifier" =~ [_]) || ("$filename" =~ [-] && "$sample_identifier" =~ [-]) ]]; then
        echo "Both Basespace sample name and input sample identifier for $filename contain matching separators..."
        SAMPLENAME_RENAMED="$sample_identifier"
      fi
      if [[ ! ("$filename" =~ [_]) || ! ("$filename" =~ [-]) ]]; then
        echo "Filename doesn't use underscore or hyphen separators, using input sample identifier as-is"
        SAMPLENAME_RENAMED="$sample_identifier"
      fi
    done

    echo "Renamed identifier: $SAMPLENAME_RENAMED"

    #Combine non-empty read files into single file without BaseSpace filename cruft
    ##FWD Read
    lane_count=0
    for fwd_read in ./dataset_*/${SAMPLENAME_RENAMED}_*R1_*.fastq.gz; do
      if [[ -s $fwd_read ]]; then
        echo "cat fwd reads: cat $fwd_read >> ~{sample_name}_R1.fastq.gz" 
        cat $fwd_read >> ~{sample_name}_R1.fastq.gz
        lane_count=$((lane_count+1))
      fi
    done
    ##REV Read
    for rev_read in ./dataset_*/${SAMPLENAME_RENAMED}_*R2_*.fastq.gz; do
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
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}