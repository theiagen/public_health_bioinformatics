version 1.0

task dnaapler_all {
  input {
    File input_fasta
    String samplename
    Int threads = 4
    Int disk_size = 100
    Int memory = 16
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dnaapler:0.8.0"
  }
  command <<< 
    set -euo pipefail

    # Check input FASTA is valid
    echo "Validating input FASTA file..."
    if ! grep -q "^>" ~{input_fasta}; then
      echo "ERROR: Input file ~{input_fasta} is not in FASTA format." >&2
      exit 1
    else
      echo "Input FASTA file is valid."
    fi

    # dnaapler version
    dnaapler --version | tee VERSION

    # Create a subdirectory for dnaapler outputs
    output_dir="dnaapler_output"
    mkdir -p "$output_dir"
    echo "Output directory created: $output_dir"

    # Run dnaapler with the 'all' subcommand
    echo "Running dnaapler..."
    dnaapler all \
      -i ~{input_fasta} \
      -o "$output_dir" \
      -p ~{samplename} \
      -t ~{threads} \
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
    cpu: threads
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 1
    preemptible: 0
  }
}