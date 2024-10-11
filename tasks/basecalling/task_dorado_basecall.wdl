version 1.0

task basecall {
  input {
    Array[File] input_files
    String dorado_model
    String kit_name
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dorado:0.8.0"  # Default Docker image
  }

  command <<< 
    set -e

    # Define the output folders
    bam_output="output/bam/"
    fastq_output="output/fastq/"
    mkdir -p $bam_output
    mkdir -p $fastq_output

    input_files_array=(~{sep=" " input_files})

    echo "### Starting basecalling process ###"
    echo "Input files: ${input_files_array[@]}"
    echo "Dorado model: ~{dorado_model}"
    echo "Kit name: ~{kit_name}"
    echo "Output BAM directory: $bam_output"
    echo "Output FASTQ directory: $fastq_output"

    # Step 1: Basecall to BAM
    echo "### Step 1: Running dorado basecaller ###"
    for file in ${input_files_array[@]}; do
      base_name=$(basename $file .pod5)
      bam_file="${bam_output}/${base_name}.bam"
      
      echo "Processing file: $file"
      echo "Output BAM file: $bam_file"

      dorado basecaller \
        /dorado_models/~{dorado_model} \
        "$file" \
        --kit-name ~{kit_name} \
        --emit-sam \
        --output-dir $bam_output \
        --verbose || { echo "ERROR: Dorado basecaller failed for $file" >&2; exit 1; }

      echo "Basecalling completed for $file"
    done

    echo "### All basecalling steps completed ###"
    ls -lh $bam_output

    # Step 2: Demultiplex with dorado demux
    echo "### Step 2: Running dorado demux ###"
    for bam_file in ${bam_output}/*.bam; do
      demuxed_output="${fastq_output}/$(basename $bam_file .bam)"
      mkdir -p $demuxed_output

      echo "Demultiplexing BAM file: $bam_file"
      echo "Demux output directory: $demuxed_output"

      dorado demux \
        --input $bam_file \
        --output-dir $demuxed_output \
        --kit-name ~{kit_name} \
        --emit-fastq \
        --verbose || { echo "ERROR: Dorado demux failed for $bam_file" >&2; exit 1; }

      echo "Demultiplexing completed for $bam_file"
    done

    echo "### All demultiplexing steps completed ###"
    ls -lh $fastq_output

    # Step 3: Combine FASTQ files for each barcode
    echo "### Step 3: Combining FASTQ files for each barcode ###"
    find ${fastq_output} -name "*.fastq" > all_fastq_files.txt
    cat all_fastq_files.txt

    barcodes=$(grep -o 'barcode[0-9]*' all_fastq_files.txt | sort | uniq)
    echo "Found barcodes: $barcodes"

    # Concatenate FASTQ files for each barcode
    for barcode in $barcodes; do
      combined_fastq="${fastq_output}/${barcode}_combined.fastq"
      echo "Combining FASTQ files for barcode: $barcode"

      : > "$combined_fastq"
      grep "$barcode" all_fastq_files.txt | while read -r fastq_file; do
        echo "Adding $fastq_file to $combined_fastq"
        cat "$fastq_file" >> "$combined_fastq"
      done
    done

    # Handle unclassified reads
    echo "### Handling unclassified reads ###"
    unclassified_fastq="${fastq_output}/unclassified_combined.fastq"
    : > "$unclassified_fastq"
    grep -L 'barcode[0-9]*' all_fastq_files.txt | while read -r fastq_file; do
      echo "Adding unclassified $fastq_file to $unclassified_fastq"
      cat "$fastq_file" >> "$unclassified_fastq"
    done

    echo "### Final output directory structure ###"
    ls -lh $fastq_output

    # Step 4: Clean up BAM files
    echo "### Step 4: Cleaning up BAM files ###"
    rm -rf $bam_output
    echo "BAM files removed"

  >>>

  output {
    Array[File] combined_fastqs = glob("output/fastq/*_combined.fastq")
  }

  runtime {
    docker: docker
    cpu: 8
    memory: "32GB"
    gpuCount: 1
    gpuType: "nvidia-tesla-t4"
    maxRetries: 3
  }
}
