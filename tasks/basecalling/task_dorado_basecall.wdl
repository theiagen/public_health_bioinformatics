version 1.0

task dorado_basecall {
  input {
    Array[File] pod5_files
    String dorado_model = "sup" # options: "fast", "hac", "sup", or explicit model name
    String kit_name

    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.9.0-cuda12.2.0"
    Int memory = 32
  }
  command <<< 
    set -euo pipefail

    dorado --version > DORADO_VERSION 2>&1
    echo "Captured Dorado version: " $(cat DORADO_VERSION)
 
    # if user provides fast, hac, or sup, pass those strings to dorado basecaller command
    if [ "~{dorado_model}" == "fast" ] || [ "~{dorado_model}" == "hac" ] || [ "~{dorado_model}" == "sup" ]; then
      dorado_model_variable=~{dorado_model}
      echo "DEBUG: dorado_model_variable is set to: $dorado_model_variable"
    # if user provides explicit model name, for example "dna_r10.4.1_e8.2_400bps_fast@v5.0.0", then provide hardcoded path to directory in container filesystem that has the dorado models
    else
      dorado_model_variable="/dorado_models/~{dorado_model}"
      echo "DEBUG: dorado_model_variable is set to: $dorado_model_variable"
    fi

    # get path to localized pod5 files
    pod5s=(~{sep=" " pod5_files})
    pod5_dir=$(dirname ${pod5s[0]})

    # Run Dorado basecaller and log output
    dorado basecaller \
      "${dorado_model_variable}" \
      ${pod5_dir}/ \
      --kit-name ~{kit_name} \
      --no-trim \
      --output-dir "output/" \
      --verbose 2> >(tee -a dorado_basecall.log >&2) || { echo "ERROR: Dorado basecaller failed"; exit 1; }

    # Log the resolved model name
    echo "DEBUG: Parsing model name from dorado log or capturing string from user input..."
    if [ "~{dorado_model}" == "fast" ] || [ "~{dorado_model}" == "hac" ] || [ "~{dorado_model}" == "sup" ]; then
      echo "DEBUG: User provided either fast, hac, or sup as input for dorado_model variable, parsing log for explicit model name now..."
      grep -m 1 'downloading' dorado_basecall.log | sed -e 's/.*downloading //' -e 's/ with.*//' | tr -d '\n' | tee DORADO_MODEL
    else
      echo "~{dorado_model}" | tee DORADO_MODEL
    fi

    echo "Basecalling completed. BAM file renamed to:" | tee -a "dorado_basecall.log"
  >>>
  output {
    Array[File] bam_files = glob("output/*.bam")
    String dorado_docker = docker
    String dorado_version = read_string("DORADO_VERSION")
    String dorado_model_used = read_string("DORADO_MODEL")
    File dorado_log = "dorado_basecall.log"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    gpuCount: 1
    gpuType: "nvidia-tesla-t4"  
    preemptible: 0
    maxRetries: 1
  }
}
