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
    for bam_file in ~{sep=" " bam_files}; do
      # Check if the BAM file exists
      if [ ! -f "$bam_file" ]; then
        echo "ERROR: BAM file $bam_file does not exist" >&2
        exit 1
      fi

      echo "Processing BAM file: $bam_file"

      # Get a clean base name for output (removing unnecessary prefixes)
      base_name=$(basename "$bam_file" .bam | sed 's/^[^_]*_//')

      echo "Demultiplexing BAM file: $bam_file"
      echo "Output directory: Current working directory (.)"

      # Run the Dorado demux command, outputting directly to the working directory
      dorado demux \
        "$bam_file" \
        --output-dir . \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --verbose > "${base_name}_demux.log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          cat "${base_name}_demux.log" >&2
          exit 1
        }

      echo "Demultiplexing completed for $bam_file. Log: ${base_name}_demux.log"

      # Rename FASTQ files to remove unwanted prefixes
      for fastq_file in ./*.fastq; do
        new_name="${base_name}_$(basename "$fastq_file" | sed 's/^[^_]*_//')"
        mv "$fastq_file" "$new_name"
        echo "Renamed $fastq_file to $new_name"
      done

      # Gzip the renamed FASTQ files
      echo "Compressing FASTQ files for $bam_file"
      gzip ./*.fastq || {
        echo "ERROR: Failed to gzip FASTQ files for $bam_file" >&2
        exit 1
      }

      echo "FASTQ files compressed for $bam_file."
      echo "Listing output files:"
      ls -lh .
    done

    echo "### Dorado demux process completed ###"
  >>>

  output {
    Array[File] fastq_files = glob("*.fastq.gz")
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16GB"
    maxRetries: 0
  }
}

