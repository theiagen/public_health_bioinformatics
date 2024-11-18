version 1.0

task cat_lanes {
  input {
    String samplename
    
    File read1_lane1
    File read1_lane2
    File? read1_lane3
    File? read1_lane4

    File? read2_lane1
    File? read2_lane2
    File? read2_lane3
    File? read2_lane4

    Int cpu = 2
    Int disk_size = 50
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.2"
    Int memory = 4
  }
  meta {
    volatile: true
  }
  command <<<
    # exit task if anything throws an error (important for proper gzip format)
    set -euo pipefail

    # move reads into single directory
    mkdir -v reads
    mv -v ~{read1_lane1} \
          ~{read2_lane1} \
          ~{read1_lane2} \
          ~{read2_lane2} \
          ~{read1_lane3} \
          ~{read2_lane3} \
          ~{read1_lane4} \
          ~{read2_lane4} \
      reads/

    # check for valid gzipped format (this task assumes FASTQ files are gzipped - they should be coming from ILMN instruments)
    gzip -t reads/*.gz

    # run concatenate script and send STDOUT/ERR to STDOUT
    # reminder: script will skip over samples that only have R1 file present
    # reminder: script REQUIRES standard illumina file endings like: _L001_R1_001.fastq.gz and _L002_R2_001.fastq.gz
    # see script here: https://github.com/theiagen/utilities/blob/main/scripts/concatenate-across-lanes.sh
    concatenate-across-lanes.sh reads/

    # ensure newly merged FASTQs are valid gzipped format
    gzip -t reads/*merged*.gz

    # determine output filenames for outputs
    mv -v reads/*_merged_R1.fastq.gz reads/~{samplename}_merged_R1.fastq.gz
    mv -v reads/*_merged_R2.fastq.gz reads/~{samplename}_merged_R2.fastq.gz
  >>>
  output {
    File read1_concatenated = "reads/~{samplename}_merged_R1.fastq.gz"
    File? read2_concatenated = "reads/~{samplename}_merged_R2.fastq.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}