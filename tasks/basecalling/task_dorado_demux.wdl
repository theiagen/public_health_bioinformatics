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

      # Run Dorado demux command
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
    done

    echo "### Listing FASTQ files after demux ###"
    ls -lh demux_output/*.fastq

    # Process and rename FASTQ files
    for fastq_file in demux_output/*.fastq; do
      if [[ "$fastq_file" == *"unclassified"* ]]; then
        final_fastq="${kit_name}_unclassified.fastq"
      else
        barcode=$(echo "$fastq_file" | sed -E 's/.*_(barcode[0-9]+)\.fastq/\1/')
        final_fastq="${kit_name}_${barcode}.fastq"
      fi

      # If the target FASTQ exists, append; otherwise, move it
      if [ -f "$final_fastq" ]; then
        echo "Appending $fastq_file to $final_fastq"
        cat "$fastq_file" >> "$final_fastq"
        rm "$fastq_file"
      else
        echo "Creating new FASTQ: $final_fastq"
        mv "$fastq_file" "$final_fastq"
      fi
    done

    echo "### Zipping all FASTQ files ###"
    for fastq in ${kit_name}_*.fastq; do
      gzip -f "$fastq"
    done

    echo "### Dorado demux process completed ###"
    ls -lh ${kit_name}_*.fastq.gz
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
