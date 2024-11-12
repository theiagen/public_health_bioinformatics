version 1.0

task basecall {
  input {
    File input_file                    # Single POD5 file for scatter processing
    String dorado_model = "sup"         # Default model to 'sup', can be overridden with full model name see docs
    String kit_name                    # Sequencing kit name
    Int disk_size = 100
    Int memory = 32
    Int cpu = 8
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -euo pipefail

    # Create a unique output directory for each scatter job
    base_name=$(basename "~{input_file}" .pod5)
    sam_output="output/sam_${base_name}/"
    mkdir -p "$sam_output"

    echo "### Starting basecalling for ~{input_file} ###"

    # Set SAM file path with unique naming based on POD5 basename
    sam_file="$sam_output/${base_name}.sam"

    echo "Processing ~{input_file}, expected output: $sam_file"

    # Run Dorado basecaller
    dorado basecaller \
      "~{dorado_model}" \
      "~{input_file}" \
      --kit-name ~{kit_name} \
      --emit-sam \
      --no-trim \
      --output-dir "$sam_output" \
      --verbose || { echo "ERROR: Dorado basecaller failed for ~{input_file}"; exit 1; }

    # Rename the generated SAM file to the unique name based on input_file
    generated_sam=$(find "$sam_output" -name "*.sam" | head -n 1)
    mv "$generated_sam" "$sam_file"

    echo "Basecalling completed for ~{input_file}. SAM file renamed to: $sam_file"
  >>>
  
  output {
    Array[File] sam_files = glob("output/sam_${base_name}/*.sam")
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

