version 1.0

task basecall {
  input {
    Array[File] input_files                # Input POD5 files for basecalling
    Boolean use_auto_model = true          # Boolean to choose automatic model selection
    String? model_speed = "sup"            # Optional model speed: "fast", "hac", or "sup"
    String? dorado_model                   # Optional: Specific model path or version
    String kit_name                        # Kit name used for sequencing
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Create output directory for SAM files
    sam_output="output/sam/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling ###"

    # Basecalling loop for each input file
    for file in ~{sep=" " input_files}; do
      base_name=$(basename "$file" .pod5)
      sam_file="$sam_output/${base_name}.sam"

      echo "Processing $file, output: $sam_file"

      # Run Dorado basecaller based on user input
      if ~{use_auto_model} && [[ -z "~{dorado_model}" ]]; then
        echo "Using automatic model selection with speed: ~{model_speed}"
        dorado basecaller \
          --model ~{model_speed} \
          "$file" \
          --kit-name ~{kit_name} \
          --emit-sam \
          --no-trim \
          --output-dir "$sam_output" \
          --verbose || { echo "ERROR: Dorado basecaller failed for $file"; exit 1; }

      else
        echo "Using specified model: ~{dorado_model}"
        dorado basecaller \
          /dorado_models/~{dorado_model} \
          "$file" \
          --kit-name ~{kit_name} \
          --emit-sam \
          --no-trim \
          --output-dir "$sam_output" \
          --verbose || { echo "ERROR: Dorado basecaller failed for $file"; exit 1; }
      fi

      echo "Basecalling completed for $file. SAM file: $sam_file"
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
