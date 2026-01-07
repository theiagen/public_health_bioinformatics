version 1.0

task trimmomatic_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/trimmomatic:0.40"
    Int trimmomatic_min_length = 75
    Int trimmomatic_window_size = 4
    Int trimmomatic_window_quality = 30
    String? trimmomatic_override_args #Note that trimming steps occur in the same order that they are given on the command line

    Boolean trimmomatic_trim_adapters = false
    File? trimmomatic_adapter_fasta
    String? trimmomatic_adapter_trim_args

    Int disk_size = 100
    Int memory = 8
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    # set default adapter trimming file
    if [[ "~{trimmomatic_trim_adapters}" == "true" ]]; then
      if [[ -n "~{trimmomatic_adapter_fasta}" ]]; then
        ADAPTER_FILE="~{trimmomatic_adapter_fasta}"
        echo "Using user supplied adapter FASTA file for adapter trimming: '$ADAPTER_FILE'"
      else
        ADAPTER_FILE="TruSeq3-PE-2.fa"
        echo "Using default file for adapter trimming: '$ADAPTER_FILE'"
      fi

      # set default adapter trimming arguments
      if [[ -n "~{trimmomatic_adapter_trim_args}" ]]; then
        ADAPTER_TRIM_ARGS="~{trimmomatic_adapter_trim_args}"
        echo "Using user supplied adapter trimming arguments: '$ADAPTER_TRIM_ARGS'"
      else
        ADAPTER_TRIM_ARGS="2:30:10"
        echo "Using default adapter trimming arguments: '$ADAPTER_TRIM_ARGS'"
      fi

      # construct adapter trimming command
      ADAPTER_TRIM_COMMAND="ILLUMINACLIP:${ADAPTER_FILE}:${ADAPTER_TRIM_ARGS}"
    else
      echo "Adapter trimming disabled"
      ADAPTER_TRIM_COMMAND=""
    fi

    # set default trimming arguments
    TRIMMOMATIC_ARGS="SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_window_quality} MINLEN:~{trimmomatic_min_length}"

    trimmomatic PE \
      ~{read1} \
      ~{read2} \
      -baseout ~{samplename}.fastq.gz \
      ${ADAPTER_TRIM_COMMAND} \
      ~{if defined(trimmomatic_override_args) then '~{trimmomatic_override_args}' else '$TRIMMOMATIC_ARGS'} \
      &> ~{samplename}.trim.stats.txt

  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    String trimmomatic_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task trimmomatic_se {
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/trimmomatic:0.40"
    Int trimmomatic_min_length = 25
    Int trimmomatic_window_size = 4
    Int trimmomatic_window_quality = 30
    String? trimmomatic_override_args #Note that trimming steps occur in the same order that they are given on the command line

    Boolean trimmomatic_trim_adapters = false
    File? trimmomatic_adapter_fasta
    String? trimmomatic_adapter_trim_args

    Int disk_size = 100
    Int memory = 8
    Int cpu = 4
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    # set default adapter trimming file
    if [[ "~{trimmomatic_trim_adapters}" == "true" ]]; then
      if [[ -n "~{trimmomatic_adapter_fasta}" ]]; then
        ADAPTER_FILE="~{trimmomatic_adapter_fasta}"
        echo "Using user supplied adapter FASTA file for adapter trimming: '$ADAPTER_FILE'"
      else
        ADAPTER_FILE="TruSeq3-SE.fa"
        echo "Using default file for adapter trimming: '$ADAPTER_FILE'"
      fi

      # set default adapter trimming arguments
      if [[ -n "~{trimmomatic_adapter_trim_args}" ]]; then
        ADAPTER_TRIM_ARGS="~{trimmomatic_adapter_trim_args}"
        echo "Using user supplied adapter trimming arguments: '$ADAPTER_TRIM_ARGS'"
      else
        ADAPTER_TRIM_ARGS="2:30:10"
        echo "Using default adapter trimming arguments: '$ADAPTER_TRIM_ARGS'"
      fi

      # construct adapter trimming command
      ADAPTER_TRIM_COMMAND="ILLUMINACLIP:${ADAPTER_FILE}:${ADAPTER_TRIM_ARGS}"
    else
      echo "Adapter trimming disabled"
      ADAPTER_TRIM_COMMAND=""
    fi

    # set default trimming arguments
    TRIMMOMATIC_ARGS="SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_window_quality} MINLEN:~{trimmomatic_min_length}"

    trimmomatic SE \
      ~{read1} \
      ~{samplename}_trimmed.fastq.gz \
      ${ADAPTER_TRIM_COMMAND} \
      ~{if defined(trimmomatic_override_args) then '~{trimmomatic_override_args}' else '$TRIMMOMATIC_ARGS'} \
      &> ~{samplename}.trim.stats.txt

  >>>
  output {
    File read1_trimmed = "~{samplename}_trimmed.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    String trimmomatic_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}