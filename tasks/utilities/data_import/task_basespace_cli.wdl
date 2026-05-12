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

    # extract resolved name and ID from the result
    resolved_name=$(echo "${result}" | awk -F',' '{print $1}')
    resolved_id=$(echo "${result}" | awk -F',' '{print $2}')
    echo "Resolved BaseSpace ${bs_app} [name, id]: ['${resolved_name}', '${resolved_id}']"

    # list datasets within the resolved run/project and substring-match against the Name column to find the dataset ID for the sample of interest
    dataset_list=$(${bs_command} list dataset "${bs_dataset_flag}=${resolved_id}" --like-type="common.fastq" --template='{{.Name}},{{.Id}}')
    echo "Looking for sample name: '~{basespace_sample_name}'"
    echo -e "Data sets in this BaseSpace ${bs_app} include:\n...\n${dataset_list}\n..." | tr ',' '\t'

    # first attempt an exact match against the Name column
    # if that fails, attempt a substring-match against the Name column. Will error if no matches or multiple matches are found in either case
    # NOTE: this approach is not perfect and could unintentionally match with other samples if the sample name provided is not specific enough
    # ex) `sample1` matches both `sample1_L001` and `sample10_L001`
    match=$(echo "${dataset_list}" | awk -F',' '$1 == "~{basespace_sample_name}"')
    if [[ -z "${match}" ]]; then
      match=$(echo "${dataset_list}" | awk -F',' '$1 ~ "~{basespace_sample_name}"')
    fi

    # exclude blank lines from the count
    match_count=$(echo "${match}" | awk 'NF' | wc -l)

    if (( match_count == 0 )); then
      echo "ERROR: no dataset matched '~{basespace_sample_name}'" >&2
      exit 1
    elif (( match_count > 1 )); then
      echo "ERROR: multiple datasets matched '~{basespace_sample_name}' in BaseSpace ${bs_app}:" >&2
      echo "${match}" | tr ',' '\t' >&2
      exit 1
    else
      # extract the exact dataset name/ID from the the matched line
      dataset_name=$(echo "${match}" | awk -F',' '{print $1}')
      dataset_id=$(echo "${match}" | awk -F',' '{print $2}')
      echo "Found '~{basespace_sample_name}' in match: ${dataset_name} (${dataset_id})"
    fi

    # download reads by dataset ID
    mkdir -p ./dataset_${dataset_id}
    ${bs_command} download dataset \
      --id="${dataset_id}" \
      --output="./dataset_${dataset_id}" \
      --extension=".fastq.gz"

    # BaseSpace will sometimes replace underscores with hyphens in the sample name
    # create a flexible pattern match to find and output the downloaded FASTQ files
    flexible_sample_name=$(echo "~{basespace_sample_name}" | sed 's/[-_]/[-_]/g')

    #Combine non-empty read files into single file without BaseSpace filename cruft
    ##FWD Read
    lane_count=0
    for fwd_read in $(find ./dataset_* -name "*${flexible_sample_name}*_R1_*.fastq.gz"); do
      if [[ -s $fwd_read ]]; then
        echo "Concatenating forward reads: ${fwd_read} >> ~{basespace_sample_name}_R1.fastq.gz"
        cat "${fwd_read}" >> "~{basespace_sample_name}_R1.fastq.gz"
        lane_count=$((lane_count+1))
      fi
    done
    ##REV Read
    for rev_read in $(find ./dataset_* -name "*${flexible_sample_name}*_R2_*.fastq.gz"); do
      if [[ -s $rev_read ]]; then
        echo "Concatenating reverse reads: ${rev_read} >> ~{basespace_sample_name}_R2.fastq.gz"
        cat "${rev_read}" >> "~{basespace_sample_name}_R2.fastq.gz"
      fi
    done
    echo "Lane Count: ${lane_count}"
  >>>
  output {
    File read1 = "~{basespace_sample_name}_R1.fastq.gz"
    File? read2 = "~{basespace_sample_name}_R2.fastq.gz"
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