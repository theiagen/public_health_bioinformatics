version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Create output directory for SAM files
    sam_output="output/sam/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling ###"

    echo "### About to list all files and directories ###"
    find /cromwell_root -type f -exec ls -lh {} \;
    echo "### Finished listing files and directories ###"


    # Basecalling loop for each input file
    for file in ~{sep=" " input_files}; do
      base_name=$(basename "$file" .pod5)
      sam_file="$sam_output/${base_name}.sam"

      echo "Processing $file, output: $sam_file"

      # Run Dorado basecaller and continue to the next file if it fails
      if dorado basecaller \
        /dorado_models/~{dorado_model} \
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
    maxRetries: 0
  }
}
