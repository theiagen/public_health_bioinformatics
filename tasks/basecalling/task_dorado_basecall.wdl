version 1.0

task basecall {
  input {
    File input_file                    # Single POD5 file for scatter processing
    String dorado_model = "sup"         # Default model to 'sup', can be overridden with full model name see docs
    String kit_name                    # Sequencing kit name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Create output directory for SAM files
    sam_output="output/sam/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling for ~{input_file} ###"

    base_name=$(basename "~{input_file}" .pod5)
    sam_file="$sam_output/${base_name}.sam"

    echo "Processing ~{input_file}, output: $sam_file"

    # Run Dorado basecaller
    dorado basecaller \
      "~{dorado_model}" \
      "~{input_file}" \
      --kit-name ~{kit_name} \
      --emit-sam \
      --no-trim \
      --output-dir "$sam_output" \
      --verbose || { echo "ERROR: Dorado basecaller failed for ~{input_file}"; exit 1; }

    echo "Basecalling completed for ~{input_file}. SAM file: $sam_file"
  >>>

  output {
    Array[File] sam_files = glob("output/sam/*.sam")
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "32GB"
    gpuCount: 1
    gpuType: "nvidia-tesla-t4"  }
}
