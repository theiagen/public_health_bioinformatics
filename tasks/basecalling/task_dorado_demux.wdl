version 1.0

task dorado_demux {
  input {
    Array[File] bam_files
    String kit_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"
  }

  command <<< 
    set -e

    # Define output folder for FASTQ files
    fastq_output="output/fastq/"
    mkdir -p "$fastq_output"

    echo "### Running Dorado demux ###"
    for bam_file in ~{sep=" " bam_files}; do
      # Logging the BAM file details
      if [ ! -f "$bam_file" ]; then
        echo "ERROR: BAM file $bam_file does not exist" >&2
        exit 1
      fi
      echo "BAM file $bam_file found. Proceeding with demux."

      base_name=$(basename "$bam_file" .bam)
      demux_output="${fastq_output}/${base_name}"
      mkdir -p "$demux_output"

      # Verify output directory creation
      if [ ! -d "$demux_output" ]; then
        echo "ERROR: Demux output directory $demux_output was not created." >&2
        exit 1
      fi

      echo "Demultiplexing BAM file: $bam_file"
      echo "Demux output directory: $demux_output"

      # Run the Dorado demux command and log output to a file
      dorado demux \
        "$bam_file" \
        --output-dir "$demux_output" \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --verbose > "$demux_output/demux.log" 2>&1 || { 
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          cat "$demux_output/demux.log" >&2
          exit 1
        }

      echo "Demultiplexing completed for $bam_file. Check $demux_output/demux.log for details."

      # List contents of the demux output directory
      echo "Contents of $demux_output:"
      ls -lh "$demux_output"
    done

    echo "### All demultiplexing steps completed ###"
    echo "Listing contents of $fastq_output:"
    ls -lh "$fastq_output"
  >>>

  output {
    Array[File] fastq_files = glob("output/fastq/**/*.fastq")
  }

  runtime {
    docker: docker
    cpu: 4
    memory: "16GB"
    maxRetries: 0
  }
}
