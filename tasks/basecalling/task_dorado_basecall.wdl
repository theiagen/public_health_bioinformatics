version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Define the top-level fastq output folder
    output_base="output/fastq/"
    mkdir -p $output_base

    input_files_array=(~{sep=" " input_files})

    echo "Input files: ${input_files_array[@]}"
    echo "Dorado model: ~{dorado_model}"

    # Run Dorado basecaller on all input files
    dorado basecaller \
      /dorado_models/~{dorado_model} \
      ${input_files_array[@]} \
      --device cuda:all \
      --emit-fastq \
      --output-dir ${output_base} > ${output_base}/basecall.log 2>&1 || { echo "Dorado basecaller failed" >&2; exit 1; }

    # Verify if any FASTQ files were created
    if ! ls ${output_base}/*.fastq 1> /dev/null 2>&1; then
      echo "Error: No FASTQ files generated" >&2
      exit 1
    fi

    # Log the final directory structure for debugging
    echo "Final output directory structure:"
    ls -lh $output_base
  >>>

  output {
    # Capture all FASTQ files directly under the output folder
    Array[File] basecalled_fastqs = glob("output/fastq/*.fastq")
    Array[File] logs = glob("output/fastq/basecall.log")
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
