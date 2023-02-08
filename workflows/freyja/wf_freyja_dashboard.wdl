version 1.0

import "../tasks/task_versioning.wdl" as versioning

workflow freyja_dashboard {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String] collection_date
    Array[String] viral_load
    String freyja_dashboard_title
    File? dashboard_intro_text
  }
  call freyja_dashboard_task {
    input:
      samplename = samplename,
      freyja_demixed = freyja_demixed,
      collection_date = collection_date,
      viral_load = viral_load,
      freyja_dashboard_title = freyja_dashboard_title,
      dashboard_intro_text = dashboard_intro_text,
  }
  call versioning.version_capture{
    input:
  }
  output {
    # Version Capture
    String freyja_dashboard_wf_version = version_capture.phvg_version
    String freyja_dashboard_wf_analysis_date = version_capture.date
    # Freyja Dashboard Visualization
    String freyja_dashboard_version = freyja_dashboard_task.freyja_dashboard_version
    File freyja_dasbhoard = freyja_dashboard_task.freyja_dasbhoard
    File freyja_demixed_aggregate = freyja_dashboard_task.freyja_demixed_aggregate
    File freyja_dashboard_metadata = freyja_dashboard_task.freyja_dashboard_metadata
  }
}

task freyja_dashboard_task {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String] collection_date
    Array[String] viral_load
    File? config
    Float? thresh
    Float? mincov
    String? headerColor
    Boolean scale_by_viral_load = false
    String freyja_dashboard_title
    File? dashboard_intro_text
    String docker = "quay.io/staphb/freyja:1.3.10"
    Int disk_size = 100
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION

  # create bash arrays
  freyja_demixed_array="~{sep=' ' freyja_demixed}"
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")
  collection_date_array=(~{sep=' ' collection_date})
  collection_date_array_len=$(echo "${#collection_date_array[@]}")
  viral_load_array=(~{sep=' ' viral_load})
  viral_load_array_len=$(echo "${#viral_load_array[@]}")

  if [ "$samplename_array_len" -ne "$collection_date_array_len" ] ||  [ "$samplename_array_len" -ne "$viral_load_array_len" ]; then
    echo "ERROR: Missing collection date or viral load value. Samplename array (length: $samplename_array_len), collection date array (length: $collection_date_array_len), and viral load array (length: $viral_load_array_len) are of unequal length." >&2
    exit 1
  else
    echo "Samplename array (length: $samplename_array_len), collection date array (length: $collection_date_array_len), and viral load array (length: $viral_load_array_len) are  of equal length." >&2.
  fi

  echo "Sample,sample_collection_datetime,viral_load" > freyja_dash_metadata.csv

    for index in ${!samplename_array[@]}; do
      samplename=${samplename_array[$index]}
      collection_date=${collection_date_array[$index]}
      viral_load=${viral_load_array[$index]}
      echo "${samplename},${collection_date},${viral_load}" >> freyja_dash_metadata.csv
    done

  # move all assemblies into single directory and aggregate files
  mkdir ./demixed_files/
  echo "mv ${freyja_demixed_array[@]} demixed_files/"
  mv ${freyja_demixed_array[@]} ./demixed_files/

  freyja aggregate \
      ./demixed_files/ \
      --output demixed_aggregate.tsv

  # Create title file
  echo "~{freyja_dashboard_title}" > dashboard-title.txt

  # Create intro text file
  if [[ ! -z "~{dashboard_intro_text}" ]]; then 
    cp "~{dashboard_intro_text}" introContent.txt
  else
    echo "SARS-CoV-2 lineage de-convolution performed by the Freyja workflow (https://github.com/andersen-lab/Freyja)." > introContent.txt
  fi

  # create freya dashboard
  echo "Running: freyja dash \
    demixed_aggregate.tsv \
    freyja_dash_metadata.csv \
    dashboard-title.txt \
    introContent.txt \
    ~{'--config ' + config} \
    ~{'--thresh ' + thresh} \
    ~{'--headerColor ' + headerColor} \
    ~{'--mincov ' + mincov} \
    ~{true='--scale_by_viral_load' false='' scale_by_viral_load} \
    --output ~{freyja_dashboard_title}.html"
  freyja dash \
    demixed_aggregate.tsv \
    freyja_dash_metadata.csv \
    dashboard-title.txt \
    introContent.txt \
    ~{'--config ' + config} \
    ~{'--thresh ' + thresh} \
    ~{'--headerColor ' + headerColor} \
    ~{'--mincov ' + mincov} \
    ~{true='--scale_by_viral_load' false='' scale_by_viral_load} \
    --output ~{freyja_dashboard_title}.html
  >>>
  output {
    String freyja_dashboard_version = read_string("FREYJA_VERSION")
    File freyja_dasbhoard = "~{freyja_dashboard_title}.html"
    File freyja_demixed_aggregate = "demixed_aggregate.tsv"
    File freyja_dashboard_metadata = "freyja_dash_metadata.csv"
  }
  runtime {
    memory: "4 GB"
    cpu: 2
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
}