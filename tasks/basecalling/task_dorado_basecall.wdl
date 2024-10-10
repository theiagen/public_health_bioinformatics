version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
    String? docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"  # Default Docker image
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
      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --output-dir ${output_base} \
        --verbose || { echo "Dorado basecaller failed for $file" >&2; exit 1; }
    done

    # Output a list of all generated FASTQ files
    echo "Listing all generated FASTQ files..."
    find ${output_base} -name "*.fastq" > all_fastq_files.txt
    cat all_fastq_files.txt

    # Identify unique barcodes from the generated files
    echo "Identifying barcodes..."
    barcodes=$(grep -o 'barcode[0-9]*' all_fastq_files.txt | sort | uniq)
    echo "Found barcodes: $barcodes"

    # Concatenate FASTQ files for each barcode
    echo "Combining FASTQ files for each barcode..."
    for barcode in $barcodes; do
      combined_fastq="${output_base}/${barcode}_combined.fastq"
      
      # Clear the target file before appending
      : > "$combined_fastq"

      grep "$barcode" all_fastq_files.txt | while read -r fastq_file; do
        cat "$fastq_file" >> "$combined_fastq"
      done
    done

    # Handle unclassified reads
    echo "Combining unclassified FASTQ files..."
    unclassified_fastq="${output_base}/unclassified_combined.fastq"
    : > "$unclassified_fastq"
    grep -L 'barcode[0-9]*' all_fastq_files.txt | while read -r fastq_file; do
      cat "$fastq_file" >> "$unclassified_fastq"
    done

    # Log final output structure
    echo "Final output directory structure:"
    ls -lh $output_base
  >>>

  output {
    Array[File] combined_fastqs = glob("output/fastq/*_combined.fastq")
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
