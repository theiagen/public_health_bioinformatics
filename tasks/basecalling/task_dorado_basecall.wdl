version 1.0

task basecall {
  input {
    Array[File] input_files
    Array[String] sample_names
    String dorado_model
    String output_prefix
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e
    mkdir -p output/fastqs
    mkdir -p output/logs

    input_files_array=(~{sep=" " input_files})
    sample_names_array=(~{sep=" " sample_names})

    echo "Input files: ${input_files_array[@]}"
    echo "Sample names: ${sample_names_array[@]}"
    echo "Dorado model: ~{dorado_model}"

    # Loop over each input file with its corresponding sample name
    for i in ${!input_files_array[@]}; do
      file=${input_files_array[i]}
      sample_name=${sample_names_array[i]}

      echo "Processing file $file with sample name $sample_name"

      # Run Dorado basecaller
      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --device cuda:all \
        --emit-fastq \
        --output-dir output/fastqs > output/logs/${sample_name}_basecaller.log 2>&1 || { echo "Dorado basecaller failed for $file" >&2; exit 1; }

      # Check for any generated FASTQ file(s) for the current sample
      generated_fastqs=(output/fastqs/*.fastq)

      if [[ ${#generated_fastqs[@]} -gt 0 ]]; then
        # Rename each FASTQ file with sample name and an index
        for idx in ${!generated_fastqs[@]}; do
          mv "${generated_fastqs[$idx]}" "output/fastqs/basecalled_${sample_name}_part${idx}.fastq"
        done
        echo "FASTQ files generated for $sample_name"
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
    Array[File] basecalled_fastqs = glob("output/fastqs/*.fastq")
    Array[File] logs = glob("output/logs/*.log")
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "32GB"
    gpuCount: 1
    gpuType: "nvidia-tesla-v100"
  }
}
