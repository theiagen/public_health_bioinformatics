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
    
    # Create a temporary directory for demux output
    mkdir -p demux_output

    for bam_file in ~{sep=" " bam_files}; do
      echo "Processing BAM file: $bam_file"

      # Run the Dorado demux command and output FASTQ files to 'demux_output/'
      dorado demux \
        "$bam_file" \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --output-dir demux_output \
        --verbose > "demux_output/$(basename $bam_file .bam)_demux.log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          exit 1
        }

      echo "Demultiplexing completed for $bam_file."
    done

    echo "### Merging FASTQ files by barcode ###"

    # Merge FASTQ files by barcode
    for fastq_file in demux_output/*.fastq; do
      barcode=$(basename "$fastq_file" | cut -d'_' -f2)  # Extract barcode from filename
      merged_fastq="~{kit_name}_${barcode}.fastq"

      # Append to existing FASTQ or create new if it doesn't exist
      if [ -f "$merged_fastq" ]; then
        echo "Appending $fastq_file to $merged_fastq"
        cat "$fastq_file" >> "$merged_fastq"
      else
        echo "Creating new FASTQ: $merged_fastq"
        mv "$fastq_file" "$merged_fastq"
      fi
    done

    echo "### Zipping FASTQ files ###"
    
    # Zip all merged FASTQ files
    for merged_fastq in ~{kit_name}_*.fastq; do
      gzip "$merged_fastq"
      echo "Zipped $merged_fastq to ${merged_fastq}.gz"
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

