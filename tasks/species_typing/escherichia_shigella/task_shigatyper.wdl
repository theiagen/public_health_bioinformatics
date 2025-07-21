version 1.0

task shigatyper {
  meta {
    description: "ShigaTyper is a quick and easy tool designed to determine Shigella serotype using Illumina (single or paired-end) or Oxford Nanopore reads with low computation requirement. https://github.com/CFSAN-Biostatistics/shigatyper"
  }
  input {
    File read1 
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/shigatyper:2.0.5"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    Boolean read1_is_ont = false
  }
  command <<<
    set -euo pipefail
    
    # get version information
    shigatyper --version | sed 's/ShigaTyper //' | tee VERSION.txt

    # if read2 DOES NOT EXIST, ASSUME SINGLE END OR ONT
    if [ -z "~{read2}" ] ; then
      INPUT_READS="--SE ~{read1}"
      # if read1_is_ont is set to TRUE, then use ONT flags
      if [ "~{read1_is_ont}" == "true" ]; then
        INPUT_READS="--SE ~{read1} --ont"
      fi
    # else read2 DOES EXIST, ASSUME PAIRED END
    else
      INPUT_READS="--R1 ~{read1} --R2 ~{read2}"
    fi
    echo "INPUT_READS set to: ${INPUT_READS}"
    echo 

    # run shigatyper. 2 output files will be ~{samplename}.tsv and ~{samplename}-hits.tsv
    echo "Running ShigaTyper..."
    shigatyper \
      ${INPUT_READS} \
      -n ~{samplename}

    # rename output TSVs to be more descriptive
    mv -v ~{samplename}.tsv ~{samplename}_shigatyper_summary.tsv

    # if 0 reads map to the reference sequences, shigatyper will not produce a hits file
    # so check for the existence of the hits file before renaming
    if [ -f ~{samplename}-hits.tsv ]; then
      echo "Shigatyper hits file exists, renaming..."
      mv -v ~{samplename}-hits.tsv ~{samplename}_shigatyper_hits.tsv
    else
      echo "Hits file does not exist, creating empty hits file..."
      touch ~{samplename}_shigatyper_hits.tsv
    fi

    # parse summary tsv for prediction, ipaB absence/presence, and notes
    cut -f 2 ~{samplename}_shigatyper_summary.tsv | tail -n 1 > shigatyper_prediction.txt
    cut -f 3 ~{samplename}_shigatyper_summary.tsv | tail -n 1 > shigatyper_ipaB_presence_absence.txt
    cut -f 4 ~{samplename}_shigatyper_summary.tsv | tail -n 1 > shigatyper_notes.txt

    # if shigatyper notes field (really the txt file) is EMPTY, write string saying it is empty to float to Terra table
    if [ "$(cat shigatyper_notes.txt)" == "" ]; then
       echo "ShigaTyper notes field was empty" > shigatyper_notes.txt
    fi

  >>>
  output {
    String shigatyper_predicted_serotype = read_string("shigatyper_prediction.txt")
    String shigatyper_ipaB_presence_absence = read_string("shigatyper_ipaB_presence_absence.txt")
    String shigatyper_notes = read_string("shigatyper_notes.txt")
    File shigatyper_hits_tsv = "~{samplename}_shigatyper_hits.tsv" # A tab-delimited detailed report file
    File shigatyper_summary_tsv = "~{samplename}_shigatyper_summary.tsv" # A tab-delimited summary report file
    String shigatyper_version = read_string("VERSION.txt")
    String shigatyper_docker = docker
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