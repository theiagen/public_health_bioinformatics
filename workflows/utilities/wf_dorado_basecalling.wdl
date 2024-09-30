version 1.0

task dorado_basecall {

  meta {
    description: "This task performs Dorado basecalling on POD5 files."
  }

  input {
    Array[File] input_files
    String dorado_model
    String output_prefix
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0" 
  }

  command <<< 
    set -e
    mkdir -p output 

    # Copy input files 
    for file in ~{sep=" " input_files}; do 
      cp "$file" ./ || { echo "Error copying $file" >&2; exit 1; }
    done 

    # List input files 
    INPUT_FILES=$(ls *.pod5)

    if [[ -z "$INPUT_FILES" ]]; then
      echo "No POD5 files found" >&2
      exit 1
    fi

    # Run Dorado basecaller using GPU 
    dorado basecaller \
      ~{dorado_model} \
      $INPUT_FILES \
      --device cuda:all \
      --emit-fastq \
      --output-dir output 

    # Rename output 
    if ls output/*.fastq 1> /dev/null 2>&1; then
      mv output/*.fastq output/~{output_prefix}.fastq
    else
      echo "Error: No FASTQ output generated" >&2
      exit 1
    fi
  >>>

  output {
    File basecalled_fastq = "output/~{output_prefix}.fastq"
  }

  runtime {
    docker: "~{docker}" 
    cpu: 8
    memory: "32GB"
    gpuCount: 1
  }
}

# To run this workflow locally, you can still provide a different Docker image:
# miniwdl run /home/fraser_combe_theiagen_com/workflows/dorado_basecalling.wdl \
#   input_files=/home/fraser_combe_theiagen_com/workflows/dna_r10.4.1_e8.2_260bps-FLO_PRO114-SQK_NBD114_96_260-4000.pod5 \
#   dorado_model=/dorado_models/dna_r10.4.1_e8.2_260bps_sup@v3.5.2 \
#   output_prefix=sample_output \
#   docker=us-docker.pkg.dev/general-theiagen/staphb/dorado:new_version \
#   --debug --verbose
