version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String output_prefix
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Define the top-level fastq output folder
    output_base="output/fastq/"
    mkdir -p $output_base

    input_files_array=(~{sep=" " input_files})

    # Initialize arrays to store paths of output files
    fastq_paths=()
    log_paths=()

    echo "Input files: ${input_files_array[@]}"
    echo "Dorado model: ~{dorado_model}"

    # Loop over each input file
    for file in ${input_files_array[@]}; do
      base_name=$(basename $file .pod5)

      # Extract barcode (assuming 'barcodeXX' is part of the base name)
      barcode=$(echo $base_name | grep -o 'barcode[0-9]\+')

      # Create a directory for the current barcode inside the fastq folder
      barcode_dir="${output_base}${barcode}"
      mkdir -p ${barcode_dir}/fastqs
      mkdir -p ${barcode_dir}/logs

      echo "Processing file $file with base name $base_name"

      # Run Dorado basecaller and store outputs in the barcode-specific directories
      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --device cuda:all \
        --emit-fastq \
        --output-dir ${barcode_dir}/fastqs > ${barcode_dir}/logs/${base_name}_basecaller.log 2>&1 || { echo "Dorado basecaller failed for $file" >&2; exit 1; }

      # Rename each FASTQ file with base name and an index
      generated_fastqs=(${barcode_dir}/fastqs/*.fastq)
      if [[ ${#generated_fastqs[@]} -gt 0 ]]; then
        for idx in ${!generated_fastqs[@]}; do
          mv "${generated_fastqs[$idx]}" "${barcode_dir}/fastqs/basecalled_${base_name}_part${idx}.fastq"
          # Store the path of each FASTQ file in the array
          fastq_paths+=("${barcode_dir}/fastqs/basecalled_${base_name}_part${idx}.fastq")
        done
        echo "FASTQ files generated for $base_name"
      else
        echo "Error: No FASTQ generated for $file" >&2
        exit 1
      fi

      # Store the path of the log file
      log_paths+=("${barcode_dir}/logs/${base_name}_basecaller.log")
    done

    # Log the final directory structure for debugging
    echo "Final output directory structure:"
    ls -lh $output_base

    # Output the list of fastq and log files
    echo "${fastq_paths[@]}" > fastq_files.txt
    echo "${log_paths[@]}" > log_files.txt
  >>>

  output {
    # Use read_lines to retrieve the paths of the generated fastqs and logs
    Array[File] basecalled_fastqs = read_lines("fastq_files.txt")
    Array[File] logs = read_lines("log_files.txt")
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
