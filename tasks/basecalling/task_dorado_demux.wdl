version 1.0

task dorado_demux {
  input {
    Array[File] bam_files
    String kit_name
    String output_file_prefix 
    Int cpu = 4 
    Int memory = 16
    Int disk_size = 100
    Boolean demux_notrim = false
    String dorado_model_used
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.9.0-cuda12.2.0"
  }
  command <<< 
    set -euo pipefail

    # Capture Dorado version
    dorado --version 2>&1 | head -n1 | tee DORADO_VERSION

    # Output the Dorado model used
    echo "~{dorado_model_used}" > DORADO_MODEL_USED

    output_file_prefix="~{output_file_prefix}"
    kit_name="~{kit_name}"

    # Start the main log file for the entire task
    exec > >(tee -a dorado_demux_output.log) 2>&1
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
        --threads ~{cpu} \
        --emit-summary \
        ~{if demux_notrim then "--no-trim" else ""} \
        --verbose > "$demux_dir/demux_${base_name}.log" 2>&1 || {
          echo "ERROR: Dorado demux failed for $bam_file. Check $demux_dir/demux_${base_name}.log for details." >&2
          exit 1
      }

      echo "Demultiplexing completed for $bam_file"
    done

    echo "### Merging FASTQ files by barcode ###"
    mkdir -p merged_output
    for demux_dir in demux_output_*; do
      for fastq_file in "$demux_dir"/*.fastq; do
        echo "Processing $fastq_file in $demux_dir"

        if [[ "$fastq_file" == *"unclassified"* ]]; then
          final_fastq="merged_output/${output_file_prefix}-unclassified.fastq"
        else
          barcode=$(echo "$fastq_file" | sed -E 's/.*_(barcode[0-9]+)\.fastq/\1/')
          final_fastq="merged_output/${output_file_prefix}-${barcode}.fastq"
        fi

        if [ -f "$final_fastq" ]; then
          cat "$fastq_file" >> "$final_fastq"
        else
          mv "$fastq_file" "$final_fastq"
        fi
      done
    done

    echo "### Compressing merged FASTQ files ###"
    pigz merged_output/*.fastq

    echo "### Dorado demux process completed successfully ###"
  >>>
  output {
    Array[File] fastq_files = glob("merged_output/*.fastq.gz")
    String dorado_docker = docker
    String dorado_version = read_string("DORADO_VERSION")
    File dorado_demux_log = "dorado_demux_output.log" 
    String dorado_model_name = read_string("DORADO_MODEL_USED")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}