version 1.0

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
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:2.0.1"
    Int disk_size = 100
    Int memory = 2
    Int cpu = 1
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
    File freyja_dashboard = "~{freyja_dashboard_title}.html"
    File freyja_demixed_aggregate = "demixed_aggregate.tsv"
    File freyja_dashboard_metadata = "freyja_dash_metadata.csv"
  }
  runtime {
    memory: memory + " GB"
    cpu: cpu
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
}
