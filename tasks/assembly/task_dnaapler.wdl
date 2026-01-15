version 1.0

task dnaapler {
  input {
    File input_fasta
    String samplename
    String dnaapler_mode = "all" # The mode of reorientation to execute (default: 'all')
    Int cpu = 4
    Int disk_size = 100
    Int memory = 16
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dnaapler:1.0.1"
  }
  command <<< 
    set -euo pipefail

    # dnaapler version
    dnaapler --version | tee VERSION

    # Create a subdirectory for dnaapler outputs
    output_dir="dnaapler_output"
    mkdir -p "$output_dir"
    echo "Output directory created: $output_dir"

    # Run dnaapler subcommand
    echo "Running dnaapler..."
    dnaapler ~{dnaapler_mode} \
      -i ~{input_fasta} \
      -o "$output_dir" \
      -p ~{samplename} \
      -t ~{cpu} \
      -f || {
        echo "ERROR: dnaapler command failed. Check logs for details." >&2
        exit 1
      }

    echo "dnaapler command completed successfully."

    # Check if output FASTA file exists
    if [[ ! -f "$output_dir"/~{samplename}_reoriented.fasta ]]; then
      echo "ERROR: Expected output file not found: $output_dir/~{samplename}_reoriented.fasta" >&2
      exit 1
    fi

    # Move the final reoriented FASTA file to the task's working directory
    echo "Moving output FASTA file to working directory..."
    mv "$output_dir"/~{samplename}_reoriented.fasta .

    echo "dnaapler task completed successfully for sample: ~{samplename}"
  >>>
  output {
    File reoriented_fasta = "~{samplename}_reoriented.fasta"
    String dnaapler_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}