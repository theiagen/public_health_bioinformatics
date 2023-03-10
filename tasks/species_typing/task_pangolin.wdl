version 1.0

task pangolin4 {
  input {
    File fasta
    String samplename
    Int min_length = 10000
    Float max_ambig = 0.5
    String docker = "staphb/pangolin:4.2-pdata-1.18.1.1"
    String? analysis_mode
    Boolean expanded_lineage=true
    Boolean skip_scorpio=false
    Boolean skip_designation_cache=false
    String? pangolin_arguments
    Int disk_size = 100
  }
  command <<<
    set -e

    # date and version capture
    date | tee DATE

    { pangolin --all-versions && usher --version; } | tr '\n' ';'  | cut -f -6 -d ';' | tee VERSION_PANGOLIN_ALL

    pangolin "~{fasta}" \
       ~{'--analysis-mode ' + analysis_mode} \
       ~{'--min-length ' + min_length} \
       ~{'--max-ambig ' + max_ambig} \
       ~{true='--expanded-lineage' false='' expanded_lineage} \
       ~{true='--skip-scorpio' false='' skip_scorpio} \
       ~{true='--skip-designation-cache' false='' skip_designation_cache} \
       --outfile "~{samplename}.pangolin_report.csv" \
       --verbose \
       ~{pangolin_arguments}


    python3 <<CODE
    import csv
    #grab output values by column header
    with open("~{samplename}.pangolin_report.csv",'r') as csv_file:
      csv_reader = list(csv.DictReader(csv_file, delimiter=","))
      for line in csv_reader:
        with open("PANGO_ASSIGNMENT_VERSION", 'wt') as assignment_version:
          pangolin_version=line["pangolin_version"]
          version=line["version"]
          assignment_version.write(f"pangolin {pangolin_version}; {version}")
        with open("PANGOLIN_LINEAGE", 'wt') as lineage:
          lineage.write(line["lineage"])
        with open("PANGOLIN_CONFLICTS", 'wt') as pangolin_conflicts:
          pangolin_conflicts.write(line["conflict"])
        with open("PANGOLIN_NOTES", 'wt') as pangolin_notes:
          pangolin_notes.write(line["note"])
        with open("EXPANDED_LINEAGE", "wt") as expanded_lineage:
          if "expanded_lineage" in line:
            expanded_lineage.write(line["expanded_lineage"])
          else:
            expanded_lineage.write("NA")
    CODE
  >>>
  output {
    String date = read_string("DATE")
    String pangolin_lineage = read_string("PANGOLIN_LINEAGE")
    String pangolin_lineage_expanded = read_string("EXPANDED_LINEAGE")
    String pangolin_conflicts = read_string("PANGOLIN_CONFLICTS")
    String pangolin_notes = read_string("PANGOLIN_NOTES")
    String pangolin_assignment_version = read_string("PANGO_ASSIGNMENT_VERSION")
    String pangolin_versions = read_string("VERSION_PANGOLIN_ALL")
    String pangolin_docker = docker
    File pango_lineage_report = "~{samplename}.pangolin_report.csv"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task pangolin_update_log {
  input {
    String samplename
    String current_lineage
    String current_pangolin_docker
    String current_pangolin_assignment_version
    String current_pangolin_versions
    String updated_lineage
    String updated_pangolin_docker
    String updated_pangolin_assignment_version
    String updated_pangolin_versions
    String? timezone
    File? lineage_log
    Int disk_size = 100
  }
  command <<<
    # set timezone for date outputs
    ~{default='' 'export TZ=' + timezone}
    DATE=$(date +"%Y-%m-%d")

    #check if lineage has been modified
    if [[ "~{current_lineage}" == "~{updated_lineage}" ]]
    then
      UPDATE_STATUS="pango lineage unchanged: ~{updated_lineage}"
    else
      UPDATE_STATUS="pango lineage modified: ~{current_lineage} -> ~{updated_lineage}"
    fi

    #if a lineage log not provided, create one with headers
    lineage_log_file="~{samplename}_pango_lineage_log.tsv"

    if [ -s "~{lineage_log}" ]
    then
      echo "Lineage log provided"

      if grep -q "previous_pangolin_assignment_version" ~{lineage_log}
      then
        mv "~{lineage_log}" ${lineage_log_file}
      else
        echo "pangolin log file provided not compatible with current PHVG version"
        exit 1
      fi
   else
     echo "Creating new lineage log file as none was provided"
     echo -e "analysis_date\tmodification_status\tprevious_lineage\tprevious_pangolin_docker\tprevious_pangolin_assignment_version\tprevious_pangolin_versions\tupdated_lineage\tupdated_pangolin_docker\tupdated_pangolin_assignment_version\tupdated_pangolin_versions" > ${lineage_log_file}
   fi

     #populate lineage log file
     echo -e "${DATE}\t${UPDATE_STATUS}\t~{current_lineage}\t~{current_pangolin_docker}\t~{current_pangolin_assignment_version}\t~{current_pangolin_versions}\t~{updated_lineage}\t~{updated_pangolin_docker}\t~{updated_pangolin_assignment_version}\t~{updated_pangolin_versions}" >> "${lineage_log_file}"

    echo "${UPDATE_STATUS} (${DATE})"  | tee PANGOLIN_UPDATE
  >>>
  output {
    String pangolin_updates = read_string("PANGOLIN_UPDATE")
    File pango_lineage_log = "~{samplename}_pango_lineage_log.tsv"
  }
  runtime {
    docker: "quay.io/theiagen/utility:1.1"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}