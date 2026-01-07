version 1.0

task freyja_plot_task {
  input {
    Array[String] samplename
    Array[File] freyja_demixed
    Array[String]? collection_date
    Boolean plot_lineages=false
    Boolean plot_time=false
    String plot_time_interval="MS"
    Int plot_day_window=14
    String freyja_plot_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:2.0.1"
    Int disk_size = 100
    Int mincov = 60
    Int memory = 2
    Int cpu = 1
  }
  command <<<
  # capture version
  freyja --version | tee FREYJA_VERSION
  
  freyja_demixed_array="~{sep=' ' freyja_demixed}"
  samplename_array=(~{sep=' ' samplename})
  samplename_array_len=$(echo "${#samplename_array[@]}")

  if ~{plot_time}; then
    # create timedate metadata sheet
    collection_date_array=(~{sep=' ' collection_date})
    collection_date_array_len=$(echo "${#collection_date_array[@]}")

    if [ "$samplename_array_len" -ne "$collection_date_array_len" ]; then
      echo "ERROR: Missing collection date. Samplename array (length: $samplename_array_len) and collection date array (length: $collection_date_array_len) are of unequal length." >&2
      exit 1
    else
      echo "Samplename array (length: $samplename_array_len) and collection date array (length: $collection_date_array_len) are of equal length." >&2.
    fi

    echo "Sample,sample_collection_datetime" > freyja_times_metadata.csv

    for index in ${!samplename_array[@]}; do
      samplename=${samplename_array[$index]}
      collection_date=${collection_date_array[$index]}
      echo "${samplename},${collection_date}" >> freyja_times_metadata.csv
    done

    plot_options="--times freyja_times_metadata.csv"

    if [ ~{plot_time_interval} == "D" ]; then
      plot_options="${plot_options} --interval D --windowsize ~{plot_day_window}"
    elif [ ~{plot_time_interval} == "MS" ]; then
      plot_options="${plot_options} --interval MS"
    else
      echo "ERROR: plot time interval value (~{plot_time_interval}) not recognized. Must be either \"D\" (days) or \"MS\" (months)" >&2
      exit 1
    fi

  fi

  # move all assemblies into single directory and aggregate files
  mkdir ./demixed_files/
  echo "mv ${freyja_demixed_array[@]} demixed_files/"
  mv ${freyja_demixed_array[@]} ./demixed_files/

  freyja aggregate \
      ./demixed_files/ \
      --output demixed_aggregate.tsv

  # create freya plot
  echo "Running: freyja plot demixed_aggregate.tsv --output ~{freyja_plot_name}.pdf ${plot_options}"
  freyja plot \
      ~{true='--lineages' false ='' plot_lineages} \
      --mincov ~{mincov} \
      demixed_aggregate.tsv \
      --output ~{freyja_plot_name}.pdf \
      ${plot_options}

  >>>
  output {
    String freyja_plot_version = read_string("FREYJA_VERSION")
    File freyja_plot = "~{freyja_plot_name}.pdf"
    File demixed_aggregate = "demixed_aggregate.tsv"
    File? freyja_plot_metadata = "freyja_times_metadata.csv"
  }
  runtime {
    memory: memory + " GB"
    cpu: cpu
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
}