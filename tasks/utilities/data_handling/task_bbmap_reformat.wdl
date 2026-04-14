version 1.0

task bbmap_reformat_interleaved{
  input {
    String samplename
    File interleaved_fastq
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bbtools:38.76"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100
  }

  command <<<
    # Pipefail set to -eu to avoid zcat pipes throwing error codes with parsing is stopped
    set -eu pipefail

    # Check if interleaved_fastq input is empty or not
    if [ "$(zcat ~{interleaved_fastq} | wc -c)" -eq 0 ]; then
      echo "ERROR: ~{interleaved_fastq} is empty after decompression."
      exit 0
    else
      echo "DEBUG: ~{interleaved_fastq} is valid."
    fi

    # Check for SRA headers. These require specific reformatting as BBTools does not recognize
    # the R1 R2 designations and will push all reads to singletons.
    echo "DEBUG Checking for SRA headers."
    FIRST_LINE=$(zcat ~{interleaved_fastq} | awk 'NR==1 {print; exit}')
    echo "DEBUG: First line: $FIRST_LINE"

    # Check for SRA specific information and structure, (SRR or ERR), SRA number, read pair designation
    # The grep is structured this way to handle the presence and absence of version numbers from the header formats. 
    # This is also replicated below in the awk command
    if echo "$FIRST_LINE" | grep -qE "@(SRR|ERR)[0-9]+.*\.(1|2)"; then
      echo "DEBUG: SRA header format detected. Changing header format from '.1 .2' to '/1 /2'"

      # Decompress, replace .1 .2 designations with /1 /2 for each header line, then compress reformatted file back into fastq.gz format
      zcat ~{interleaved_fastq} | \
        awk '/^@(SRR|ERR)[0-9]+.*\.1 / {sub(/\.1 /, "/1 ")} 
              /^@(SRR|ERR)[0-9]+.*\.2 / {sub(/\.2 /, "/2 ")} 
              {print}' | \
        gzip > ~{samplename}_reformatted_sra.fastq.gz
    else
      echo "DEBUG: No SRA header format detected. Preserving format."
    fi

    # Check if the fastq has been reformated and assign INPUT_FASTQ accordingly
    if [ -s  ~{samplename}_reformatted_sra.fastq.gz ]; then
      INPUT_FASTQ="~{samplename}_reformatted_sra.fastq.gz"
    else
      INPUT_FASTQ="~{interleaved_fastq}"
    fi
    
    echo "DEBUG: INPUT_FASTQ is $INPUT_FASTQ"
    
    # Capture exit code for checking interleaved status
    # Removing e from pipefail in order to check for interleaved status as 
    # reformat.sh will return exit code 1 if not properly interleaved
    set +e
    reformat.sh in=$INPUT_FASTQ out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t
    reformat_exit_code=$?
    echo "DEBUG: Initial reformat exit code $reformat_exit_code"
    set -e

    # Check for a non 0 exit code, meaning reads need to be repaired due to mismatched pairs
    if [ $reformat_exit_code -ne 0 ]; then 
      # Run repair.sh and reformat on corrected reads
      echo "DEBUG: Names do not appear to be correctly paired in the interleaved FASTQ file. Running repair.sh"

      # Ensure repair.sh is run correctly
      if ! repair.sh in=$INPUT_FASTQ out=repaired.fastq outs=singletons.fastq repair=t overwrite=t; then
        echo "ERROR: repair.sh has failed to correct read pairs" >&2
        exit 0
      fi

      # Check for needed repair.sh outputs prior to running reformat.sh
      if [ ! -s repaired.fastq ]; then
        echo "ERROR: repair.sh produced empty output" >&2
        exit 0
      fi

      # Run reformat.sh on corrected reads. Set overwrite to true to over write any created outputs from initial run.
      echo "DEBUG: repair.sh complete, running reformat.sh to deinterleave " 
      reformat.sh in=repaired.fastq out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t \
        overwrite=t
    fi
    
    echo "DEBUG: reformat.sh complete, compressing deinterleaved FASTQs"
    # GZIP deinterleaved FASTQ files with additional error handling for missing and empty reformat output.
    for fastq_file in "~{samplename}_deinterleaved_R1.fastq" "~{samplename}_deinterleaved_R2.fastq"; do
      if [ ! -s "$fastq_file" ]; then
        echo "ERROR: Expected output file '$fastq_file' is missing or empty" >&2
        exit 0
      fi
      gzip "$fastq_file"
    done
  >>>

  output {
    File? deinterleaved_fastq_R1 = "~{samplename}_deinterleaved_R1.fastq.gz"
    File? deinterleaved_fastq_R2 = "~{samplename}_deinterleaved_R2.fastq.gz"
    String bbmap_reformat_docker = docker
  }

  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}