version 1.0

task basecall {
  input {
    File input_file                    # Single POD5 file for scatter processing
    String dorado_model = "sup"         # Default model to 'sup', can be overridden with full model name see docs
    String kit_name                    # Sequencing kit name
    Int disk_size = 100
    Int memory = 32
    Int cpu = 8
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.3"
  }
  command <<< 
    set -euo pipefail

    # Capture Dorado version and log it
    dorado --version > DORADO_VERSION 2>&1
    echo "Captured Dorado version:" $(cat DORADO_VERSION)
 
    # if user provides fast, hac, or sup, pass those strings to dorado basecaller command
    if [ "~{dorado_model}" == "fast" ] || [ "~{dorado_model}" == "hac" ] || [ "~{dorado_model}" == "sup" ]; then
      dorado_model_variable=~{dorado_model}
      echo "DEBUG: dorado_model_variable is set to: $dorado_model_variable"
    # if user provides explicit model name, for example "dna_r10.4.1_e8.2_400bps_fast@v5.0.0", then provide hardcoded path to directory in container filesystem that has the dorado models
    else
      dorado_model_variable="/dorado_models/~{dorado_model}"
      echo "DEBUG: dorado_model_variable is set to: $dorado_model_variable"
    fi

    # Create a unique output directory for each scatter job
    base_name=$(basename "~{input_file}" .pod5)
    sam_output="output/sam_${base_name}/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling for ~{input_file} ###" | tee -a dorado_basecall.log

    # Set SAM file path with unique naming based on POD5 basename
    sam_file="$sam_output/${base_name}.sam"

    echo "Processing ~{input_file}, expected output: $sam_file" | tee -a dorado_basecall.log

    # Run Dorado basecaller and log output
    # This part "2> >(tee -a log.txt >&2)" is used to redirect STDERR to the screen AND to append the STDERR to log.txt file. 
    # Useful for troubleshooting in Terra and for parsing for important information.
    dorado basecaller \
      "${dorado_model_variable}" \
      "~{input_file}" \
      --kit-name ~{kit_name} \
      --emit-sam \
      --no-trim \
      --output-dir "$sam_output" \
      --verbose 2> >(tee -a dorado_basecall.log >&2) || { echo "ERROR: Dorado basecaller failed for ~{input_file}"; exit 1; }

    # Log the resolved model name
    echo "DEBUG: Parsing model name from dorado log or capturing string from user input..."
    if [ "~{dorado_model}" == "fast" ] || [ "~{dorado_model}" == "hac" ] || [ "~{dorado_model}" == "sup" ]; then
      echo "DEBUG: User provided either fast, hac, or sup as input for dorado_model variable, parsing log for explicit model name now..."
      grep -m 1 'downloading' dorado_basecall.log | sed -e 's/.*downloading //' -e 's/ with.*//' | tr -d '\n' | tee DORADO_MODEL

    # (else) if user provides explicit model name, just output that string, no parsing involved
    else
      echo "~{dorado_model}" | tee DORADO_MODEL
    fi

    # Rename the generated SAM file to the unique name based on input_file
    generated_sam=$(find "$sam_output" -name "*.sam" | head -n 1)
    mv "$generated_sam" "$sam_file"

    echo "Basecalling completed for ~{input_file}. SAM file renamed to: $sam_file" | tee -a "dorado_basecall.log"
  >>>
  output {
    Array[File] sam_files = glob("output/sam_*/*.sam")
    String dorado_docker = docker
    String dorado_version = read_string("DORADO_VERSION")
    String dorado_model_used = read_string("DORADO_MODEL")
    File dorado_log = "dorado_basecall.log"
  }
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    gpuCount: 1
    gpuType: "nvidia-tesla-t4"  
    preemptible: 0
    maxRetries: 1
  }
}
