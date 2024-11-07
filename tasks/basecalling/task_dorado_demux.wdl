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

    echo "### Starting Dorado demux ###"
    date
    echo "Input BAM files:"
    for bam_file in ~{sep=" " bam_files}; do echo "$bam_file"; done

    mkdir -p demux_output

    for bam_file in ~{sep=" " bam_files}; do
      echo "Processing BAM file: $bam_file"
      echo "Output directory: demux_output"
      echo "Running dorado demux for kit: ~{kit_name}"

      # Run Dorado demux command
      dorado demux \
        "$bam_file" \
        --output-dir demux_output \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --emit-summary \
        --verbose > "demux_output/$(basename "$bam_file").log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file" >&2
          cat "demux_output/$(basename "$bam_file").log" >&2
          exit 1
      }

      echo "Demultiplexing completed for $bam_file"
    done

    echo "### Listing FASTQ files after demux ###"
    ls -lh demux_output/*.fastq || echo "No FASTQ files found."

    # Debugging the file naming logic
    echo "### Renaming FASTQ files ###"
    for fastq_file in demux_output/*.fastq; do
      echo "Processing $fastq_file"
      
      if [[ "$fastq_file" == *"unclassified"* ]]; then
        final_fastq="~{fastq_file_name}-unclassified.fastq"
      else
        barcode=$(echo "$fastq_file" | sed -E 's/.*_(barcode[0-9]+)\.fastq/\1/')
        final_fastq="~{fastq_file_name}-${barcode}.fastq"
      fi

      echo "Renaming $fastq_file to $final_fastq"

      if [ -f "$final_fastq" ]; then
        echo "Appending to existing $final_fastq"
        cat "$fastq_file" >> "$final_fastq"
        rm "$fastq_file"
      else
        echo "Creating new FASTQ file: $final_fastq"
        mv "$fastq_file" "$final_fastq"
      fi
    done

    echo "### Zipping all FASTQ files ###"
    for fastq in ~{fastq_file_name}-*.fastq; do
      echo "Zipping $fastq"
      gzip -f "$fastq"
    done

    echo "### Final FASTQ Files ###"
    ls -lh ~{fastq_file_name}-*.fastq.gz || echo "No gzipped FASTQ files found."

    echo "### Dorado demux process completed successfully ###"
    date
  >>>

  output {
    Array[File] fastq_files = glob("~{fastq_file_name}-*.fastq.gz")
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
