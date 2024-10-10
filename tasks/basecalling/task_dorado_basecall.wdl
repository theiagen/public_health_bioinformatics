task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Define the output folder
    output_base="output/fastq/"
    mkdir -p $output_base

    input_files_array=(~{sep=" " input_files})

    echo "Input files: ${input_files_array[@]}"
    echo "Dorado model: ~{dorado_model}"
    echo "Kit name: ~{kit_name}"

    # Run Dorado basecaller
    for file in ${input_files_array[@]}; do
      base_name=$(basename $file .pod5)
      log_file="${output_base}/${base_name}_basecall.log"

      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --output-dir ${output_base} \
        --verbose > $log_file 2>&1 || { echo "Dorado basecaller failed for $file" >&2; exit 1; }
    done

    # Concatenate all FASTQ files for each barcode
    echo "Combining FASTQ files for each barcode..."
    find ${output_base} -name "*.fastq" | while read -r fastq_file; do
      barcode_name=$(basename "$fastq_file" | cut -d'_' -f2)
      combined_fastq="${output_base}/${barcode_name}_combined.fastq"
      cat "$fastq_file" >> "$combined_fastq"
    done

    # Log final output structure
    echo "Final output directory structure:"
    ls -lh $output_base
  >>>

  output {
    Array[File] combined_fastqs = glob("output/fastq/*_combined.fastq")
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
