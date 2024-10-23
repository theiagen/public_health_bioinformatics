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

    # Debugging: List all files in the current directory structure
    echo "### Listing all files and directories before basecalling ###"
    find /cromwell_root -type f -exec ls -lh {} \;

    for file in ~{sep=" " input_files}; do
      base_name=$(basename "$file" .pod5)
      sam_file="$sam_output/${base_name}.sam"

      echo "Processing $file, output: $sam_file"

      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --kit-name ~{kit_name} \
        --emit-sam \
        --no-trim \
        --output-dir "$sam_output" \
        --verbose || { echo "ERROR: Dorado basecaller failed for $file"; exit 1; }

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
    maxRetries: 0
  }
}
