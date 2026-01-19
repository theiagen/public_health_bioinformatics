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
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    # Capture exit code for checking interleaved status
    # Removing e from pipefail in order to check for interleaved status as 
    # reformat.sh will return exit code 1 if not properly interleaved
    set +e 
    reformat.sh in=~{interleaved_fastq} out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t
    reformat_exit_code=$?
    echo "DEBUG: Initial reformat exit code $reformat_exit_code"
    # Reset pipefail after initial reformat.sh run
    set -e

    # Check for a non 0 exit code, meaning reads need to be repaired due to mismatched pairs
    if [ $reformat_exit_code -ne 0 ]; then 
      # Run repair.sh and reformat on corrected reads
      echo "DEBUG: Names do not appear to be correctly paired in the interleaved FASTQ file. Running repair.sh"
      repair.sh in=~{interleaved_fastq} out=repaired.fastq repair=t overwrite=t
      echo "DEBUG: repair.sh complete, running reformat.sh to deinterleave " 
      reformat.sh in=repaired.fastq out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t \
        overwrite=t
    fi
    
    echo "DEBUG: reformat.sh complete, compressing deinterleaved FASTQs"
    # GZIP deinterleaved FASTQ files
    gzip *_deinterleaved_R*.fastq

  >>>

  output {
    File deinterleaved_fastq_R1 = "~{samplename}_deinterleaved_R1.fastq.gz"
    File deinterleaved_fastq_R2 = "~{samplename}_deinterleaved_R2.fastq.gz"
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