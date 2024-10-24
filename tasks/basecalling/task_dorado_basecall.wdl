version 1.0

task basecall {
  input {
    Array[File] input_files
    String? dorado_model               # Optional: Manual model input
    Boolean use_auto_model = true      # Use automatic model selection if true
    String model_accuracy = "sup"      # Default to 'sup' (most accurate model)
    String kit_name                    # Sequencing kit name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Create output directory for SAM files
    sam_output="output/sam/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling ###"

    # Determine which model to use (auto or manual)
    model_to_use=$(if [[ ~{use_auto_model} == "true" ]]; then echo ~{model_accuracy}; else echo ~{dorado_model}; fi)

    # Loop through input files and basecall
    for file in ~{sep=" " input_files}; do
      base_name=$(basename "$file" .pod5)
      sam_file="$sam_output/${base_name}.sam"

      echo "Processing $file, output: $sam_file"

      # Run Dorado basecaller
      if dorado basecaller \
        "$model_to_use" \
        "$file" \
        --kit-name ~{kit_name} \
        --emit-sam \
        --no-trim \
        --output-dir "$sam_output" \
        --verbose; then
          echo "Basecalling completed successfully for $file. SAM file: $sam_file"
      else
          echo "ERROR: Dorado basecaller failed for $file. Moving on to the next file."
      fi
    done

    echo "Basecalling steps completed."
  >>>

  output {
    Array[File] sam_files = glob("output/sam/*.sam")
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "32GB"
    gpuCount: 1
    gpuType: "nvidia-tesla-t4"
    maxRetries: 3
  }
}
