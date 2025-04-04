version 1.0

task dorado_demux {
  input {
    Array[File] bam_files
    String kit_name
    String output_file_prefix

    Boolean demux_no_trim = false
    String dorado_model_used

    Int cpu = 4 
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.9.0-cuda12.2.0"    
    Int memory = 16
  }
  command <<< 
    set -euo pipefail

    dorado --version 2>&1 | head -n1 | tee DORADO_VERSION
    echo "~{dorado_model_used}" > DORADO_MODEL_USED

    echo "DEBUG: moving bam files"
    mkdir input_bams
    for bam_file in ~{sep=" " bam_files}; do 
      echo "$bam_file" 
      mv $bam_file input_bams/
    done

    # Run Dorado demux command
    dorado demux \
      input_bams/ \
      --output-dir demuxed_fastqs \
      --kit-name ~{kit_name} \
      --emit-fastq \
      --threads ~{cpu} \
      --emit-summary \
      ~{if demux_no_trim then "--no-trim" else ""} \
      --verbose 2>&1

    echo "DEBUG: demultiplexing should've finished"

    mkdir renamed_fastqs
    for file in demuxed_fastqs/*.fastq; do
      filename=$(basename $file | rev | cut -d"_" -f1 | rev)
      mv $file renamed_fastqs/~{output_file_prefix}-${filename}
    done

    echo "### Compressing merged FASTQ files ###"
    pigz renamed_fastqs/*.fastq

    echo "### Dorado demux process completed successfully ###"
  >>>
  output {
    Array[File] fastq_files = glob("renamed_fastqs/*.fastq.gz")
    String dorado_docker = docker
    String dorado_version = read_string("DORADO_VERSION")
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