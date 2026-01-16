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
    # TODO: Need to set back to -euo and sort out checking the interleaved status. reformat.sh throws an error when not
    # properly paired so this needs to be supressed somehow. 
    set -uo pipefail

    if reformat.sh in=~{interleaved_fastq} out=stdout.fq verifypaired=t; then
      echo "DEBUG: FASTQ file appears to be correctly paired."
      reformat.sh in=~{interleaved_fastq} out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t 
    else
      echo "DEBUG: Names do not appear to be correctly paired in the interleaved FASTQ file. Running repair.sh"
      repair.sh in=~{interleaved_fastq} out=repaired.fastq repair=t 
      reformat.sh in=repaired.fastq out=~{samplename}_deinterleaved_R1.fastq \
        out2=~{samplename}_deinterleaved_R2.fastq \
        verifypaired=t
    fi
  >>>

  output {
    File deinterleaved_fastq_R1 = "~{samplename}_deinterleaved_R1.fastq"
    File deinterleaved_fastq_R2 = "~{samplename}_deinterleaved_R2.fastq"
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