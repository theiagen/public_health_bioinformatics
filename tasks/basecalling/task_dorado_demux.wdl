version 1.0

task dorado_demux {
  input {
    Array[File] bam_files
    String kit_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e
    echo "### Running Dorado demux ###"
    
    mkdir -p demux_output

    for bam_file in ~{sep=" " bam_files}; do
      echo "Processing BAM file: $bam_file"

      # Run the Dorado demux command
      dorado demux \
        "$bam_file" \
        --output-dir demux_output \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --verbose > "demux_output/$(basename "$bam_file").log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          cat "demux_output/$(basename "$bam_file").log" >&2
          exit 1
      }

      echo "Demultiplexing completed for $bam_file."

      # Process the output FASTQ files
      for fastq_file in demux_output/*.fastq; do
        barcode=$(echo "$fastq_file" | sed -E 's/.*_(barcode[0-9]+)\.fastq/\1/')

        # Determine the final output file name
        final_fastq="${kit_name}_${barcode}.fastq"

        if [ -f "$final_fastq" ]; then
          echo "Appending $fastq_file to $final_fastq"
          cat "$fastq_file" >> "$final_fastq"
        else
          echo "Creating new FASTQ: $final_fastq"
          mv "$fastq_file" "$final_fastq"
        fi
      done
    done

    echo "### Zipping all FASTQ files ###"
    for fastq in ${kit_name}_*.fastq; do
      gzip -f "$fastq"
    done

    echo "### Dorado demux process completed ###"
    ls -lh ~{kit_name}_*.fastq.gz
  >>>

  output {
    Array[File] fastq_files = glob("~{kit_name}_*.fastq.gz")
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16GB"
    maxRetries: 0
  }
}
