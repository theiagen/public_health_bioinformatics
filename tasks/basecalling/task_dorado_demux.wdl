version 1.0

task dorado_demux {
  input {
    Array[File] bam_files
    String kit_name
    String fastq_file_name
    Int cpu = 4 
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -euo pipefail

    fastq_file_name="~{fastq_file_name}"
    kit_name="~{kit_name}"

    echo "### Starting Dorado demux ###"
    date
    echo "Input BAM files:"
    for bam_file in ~{sep=" " bam_files}; do echo "$bam_file"; done

    # Process each BAM file in a unique output directory
    for bam_file in ~{sep=" " bam_files}; do
      base_name=$(basename "$bam_file" .bam)
      demux_dir="demux_output_${base_name}"
      mkdir -p "$demux_dir"

      echo "Processing BAM file: $bam_file into directory $demux_dir"
      echo "Running dorado demux for kit: $kit_name"

      # Run Dorado demux command
      dorado demux \
        "$bam_file" \
        --output-dir "$demux_dir" \
        --kit-name "$kit_name" \
        --emit-fastq \
        --emit-summary \
        --verbose > "$demux_dir/demux_${base_name}.log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          cat "$demux_dir/demux_${base_name}.log" >&2
          exit 1
      }

      echo "Demultiplexing completed for $bam_file"
    done

    # Debugging merged output and merging logic
    echo "### Merging FASTQ files by barcode ###"
    mkdir -p merged_output
    for demux_dir in demux_output_*; do
      for fastq_file in "$demux_dir"/*.fastq; do
        echo "Processing $fastq_file in $demux_dir"

        # Check if the file is "unclassified" or has a specific barcode
        if [[ "$fastq_file" == *"unclassified"* ]]; then
          final_fastq="merged_output/${fastq_file_name}-unclassified.fastq"
        else
          barcode=$(echo "$fastq_file" | sed -E 's/.*_(barcode[0-9]+)\.fastq/\1/')
          final_fastq="merged_output/${fastq_file_name}-${barcode}.fastq"
        fi

        # Verify if the merged file already exists
        if [ -f "$final_fastq" ]; then
          echo "Appending to existing $final_fastq"
          cat "$fastq_file" >> "$final_fastq"
        else
          echo "Creating new FASTQ file: $final_fastq"
          mv "$fastq_file" "$final_fastq"
        fi
      done
    done

    # Compress the final merged FASTQ files
    echo "### Compressing merged FASTQ files ###"
    for final_fastq in merged_output/*.fastq; do
      echo "Compressing $final_fastq"
      gzip -f "$final_fastq"
    done

    echo "### Final Merged FASTQ Files ###"
    ls -lh merged_output/*.fastq.gz || echo "No gzipped FASTQ files found."

    echo "### Dorado demux process completed successfully ###"
    date
  >>>

  output {
    Array[File] fastq_files = glob("merged_output/*.fastq.gz")
  }
  
  runtime {
    docker: docker
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 1
  }
}
