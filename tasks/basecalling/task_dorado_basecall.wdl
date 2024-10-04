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
    mkdir -p output/logs

    input_files_array=(~{sep=" " input_files})

    echo "Input files: ${input_files_array[@]}"
    echo "Dorado model: ~{dorado_model}"

    # Loop over each input file
    for file in ${input_files_array[@]}; do
      base_name=$(basename $file .pod5)

      # Extract barcode (assuming 'barcodeXX' is part of the base name)
      barcode=$(echo $base_name | grep -o 'barcode[0-9]\+')

      # Create a directory for the current barcode
      mkdir -p output/fastqs/${barcode}

      echo "Processing file $file with base name $base_name"

      # Run Dorado basecaller
      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --device cuda:all \
        --emit-fastq \
        --output-dir output/fastqs/${barcode} > output/logs/${base_name}_basecaller.log 2>&1 || { echo "Dorado basecaller failed for $file" >&2; exit 1; }

      # Check for any generated FASTQ file(s) for the current file
      generated_fastqs=(output/fastqs/${barcode}/*.fastq)

      if [[ ${#generated_fastqs[@]} -gt 0 ]]; then
        # Rename each FASTQ file with base name and an index
        for idx in ${!generated_fastqs[@]}; do
          mv "${generated_fastqs[$idx]}" "output/fastqs/${barcode}/basecalled_${base_name}_part${idx}.fastq"
        done
        echo "FASTQ files generated for $base_name"
      else
        echo "Error: No FASTQ generated for $file" >&2
        exit 1
      fi
    done

    # Log the final directory structure for debugging
    echo "Final output directory structure:"
    ls -lh output/fastqs/
  >>>

  output {
    Array[File] basecalled_fastqs = glob("output/fastqs/*/*.fastq")
    Array[File] logs = glob("output/logs/*.log")
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
