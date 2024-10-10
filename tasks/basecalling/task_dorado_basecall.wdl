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
    for file in ${input_files_array[@]}; do
      base_name=$(basename $file .pod5)
      log_file="${output_base}/${base_name}_basecall.log"

      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --device cuda:all \
        --emit-fastq \
        --output-dir ${output_base} > $log_file 2>&1 || { echo "Dorado basecaller failed for $file" >&2; exit 1; }
    done

    # Log the final directory structure for debugging
    echo "Final output directory structure:"
    ls -lh $output_base
  >>>

  output {
    # Output all the FASTQ files and log files
    Array[File] basecalled_fastqs = glob("output/fastq/**/*.fastq")
    Array[File] logs = glob("output/fastq/*.log")
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
