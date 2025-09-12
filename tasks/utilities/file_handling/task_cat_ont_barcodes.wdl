version 1.0

task cat_ont_barcodes {
  input {
    String input_bucket_path # include trailing slash if directory
    String output_bucket_path # include trailing slash if directory
    String file_extension = ".fastq.gz" # default extension for ONT barcodes
    File? barcode_renaming_file # file containing barcode renaming mappings
    
    Int cpu = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/ont-barcodes:0.0.2"
    Int memory = 4
  }
  meta {
    volatile: true
  }
  command <<<
    set -euo pipefail

    echo "Running concatenate-barcodes.py to concatenate ONT barcodes with extension ~{file_extension} from ~{input_bucket_path} to ~{output_bucket_path}"

    concatenate-barcodes.py \
      ~{input_bucket_path} \
      ~{output_bucket_path} \
      --file_extension ~{file_extension} \
      ~{"--map_file " + barcode_renaming_file} \
      --gcp --verbose
    
    echo "Concatenation script completed."

  >>>
  output {
    File concatenation_log = "concatenate_barcodes.log"
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