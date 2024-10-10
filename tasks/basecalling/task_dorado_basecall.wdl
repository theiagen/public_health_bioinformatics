version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Define the top-level FASTQ output folder
    output_base="/cromwell_root/output/fastq/"
    mkdir -p ${output_base}

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

    # Verify that FASTQ files are generated in expected barcode subdirectories
    if ls ${output_base}/*/*.fastq 1> /dev/null 2>&1; then
      echo "FASTQ files successfully generated."
    else
      echo "Error: No FASTQ files found in ${output_base}" >&2
      exit 1
    fi

    # Log the final directory structure for debugging
    echo "Final output directory structure:"
    ls -lhR ${output_base}
  >>>

  output {
    # Collect all the FASTQ files and log files from the output folder
    Array[File] basecalled_fastqs = glob("output/fastq/*/*.fastq")
    File log = glob("output/fastq/basecall.log")[0]
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
