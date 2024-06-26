version 1.0

task stxtyper {
  input {
    File assembly
    String samplename
    String docker = "kapsakcj/stxtyper:78754d7"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 2
  }
  command <<<
    # capture date
    date | tee DATE

    # capture version info
    stxtyper --version | tee VERSION.txt

    echo "DEBUG: running StxTyper now..."
    # run StxTyper on assembly; may need to add/remove options in the future if they change
    stxtyper \
      --nucleotide ~{assembly} \
      --name ~{samplename} \
      --output ~{samplename}_stxtyper.tsv \
      --log ~{samplename}_stxtyper.log

    # parse output TSV
    echo "DEBUG: Parsing StxTyper output TSV..."

    # check for output file with only 1 line (meaning no hits found); exit cleanly if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -eq 1 ]; then
      echo "No hits found by StxTyper" > stxtyper_hits.txt
      echo "0" > stxtyper_num_hits.txt
      echo "DEBUG: No hits found in StxTyper output TSV. Exiting task with exit code 0 now."
      
      # create empty output files
      touch stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stxA_complete_hits.txt stxtyper_stxB_complete_hits.txt stxtyper_stxA_partial_hits.txt stxtyper_stxB_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stxtyper_A_partial_B_complete.txt stx_novel_hits.txt
      # put "none" into all of them so task does not fail
      echo "None" | tee stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stxA_complete_hits.txt stxtyper_stxB_complete_hits.txt stxtyper_stxA_partial_hits.txt stxtyper_stxB_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stxtyper_A_partial_B_complete.txt stx_novel_hits.txt
      exit 0
    fi
    
    # check for output file with more than 1 line (meaning hits found); count lines & parse output TSV if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -gt 1 ]; then
      echo "Hits found by StxTyper. Counting lines & parsing output TSV now..."
      # count number of lines in output TSV (excluding header)
      wc -l < ~{samplename}_stxtyper.tsv | awk '{print $1-1}' > stxtyper_num_hits.txt
      # remove header line
      sed '1d' ~{samplename}_stxtyper.tsv > ~{samplename}_stxtyper_noheader.tsv

      ### parse output TSV ###
      # complete operons
      echo "DEBUG: Parsing complete operons..."
      awk -F'\t' -v OFS=, '$4 == "COMPLETE" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stxtyper_complete_operons.txt
      # if grep for COMPLETE fails, write "None" to file for output string
      if [[ "$(grep --silent 'COMPLETE' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]]; then
        echo "None" > stxtyper_complete_operons.txt
      fi

      # complete_novel operons
      echo "DEBUG: Parsing complete novel hits..."
      awk -F'\t' -v OFS=, '$4 == "COMPLETE_NOVEL" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stx_novel_hits.txt
      # if grep for COMPLETE_NOVEL fails, write "None" to file for output string
      if [ "$(grep --silent 'COMPLETE_NOVEL' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stx_novel_hits.txt
      fi

      # stxA complete hits
      echo "DEBUG: Parsing stxA complete hits..."
      # explanation: if "A_reference_suptype" col contains stxA AND %ID >= 98.00 AND %COV >= 100.00, print the "A_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$10 ~ "stxA.*" && $11 >= 98.00 && $12 >= 100.00 {print $10}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxA_complete_hits.txt
      # if no stxA complete hits found, write "None" to file for output string
      if [ ! -s stxtyper_stxA_complete_hits.txt ]; then
        echo "None" > stxtyper_stxA_complete_hits.txt
      fi

      # stxB complete hits
      echo "DEBUG: Parsing stxB complete hits..."
      # explanation: if "B_reference_suptype" col contains stxB AND %ID >= 98.00 AND %COV >= 100.00, print the "B_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$14 ~ "stxB.*" && $15 >= 98.00 && $16 >= 100.00 {print $14}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxB_complete_hits.txt
      # if no stxB complete hits found, write "None" to file for output string
      if [ ! -s stxtyper_stxB_complete_hits.txt ]; then
        echo "None" > stxtyper_stxA_complete_hits.txt
      fi

      # partial hits (to any gene in stx operon)
      echo "DEBUG: Parsing stxtyper partial hits"
      # explanation: if "operon" column contains "PARTIAL" (either PARTIAL or PARTIAL_CONTIG_END possible); print either "stx1" or "stx2" or "stx1,stx2"
      awk -F'\t' -v OFS=, '$4 ~ "PARTIAL.*" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_partial_hits.txt
      # if no stxA complete hits found, write "None" to file for output string
      if [ ! -s stxtyper_partial_hits.txt ]; then
        echo "None" > stxtyper_partial_hits.txt
      fi

      # stxA partial hits
      echo "DEBUG: Parsing stxA partial hits..."
      # explanation: if "A_reference_suptype" col contains stxA AND  %COV < 100.00, print the "A_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$10 ~ "stxA.*" && $12 < 100.00 {print $10}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxA_partial_hits.txt
      # if no stxA partial hits found, write "None" to file for output string
      ## TODO FIX THIS CONDITIONAL, IT DOESN'T WORK SINCE SOME BLANK FILES ARE BIGGER THAN 0 BYTES
      if [ ! -s stxtyper_stxA_partial_hits.txt ]; then
        echo "None" > stxtyper_stxA_partial_hits.txt
      fi

      # stxB partial hits
      echo "DEBUG: Parsing stxB partial hits..."
      # explanation: if "B_reference_suptype" col contains stxB AND  %COV < 100.00, print the "A_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$14 ~ "stxB.*" && $16 < 100.00 {print $14}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxB_partial_hits.txt
      # if no stxB partial hits found, write "None" to file for output string
      ## TODO FIX THIS CONDITIONAL, IT DOESN'T WORK
      if [ ! -s stxtyper_stxB_partial_hits.txt ]; then
        echo "None" > stxtyper_stxB_partial_hits.txt
      fi

      # frameshifts or internal stop codons in stx genes
      echo "DEBUG: Parsing stx frameshifts or internal stop codons..."
      # explanation: if operon column contains "FRAME_SHIFT" or "INTERNAL_STOP", print the "operon" in a sorted/unique list
      awk -F'\t' -v OFS=, '$4 == "FRAME_SHIFT" || $4 == "INTERNAL_STOP" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      # if no frameshifts or internal stop codons found, write "None" to file for output string
      if [ "$(grep --silent -E 'FRAME_SHIFT|INTERNAL_STOP' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      fi

      # temporarily making a bunch of empty files to test parsing
      touch stxtyper_A_partial_B_complete.txt

    fi

    echo "DEBUG: Finished parsing StxTyper output TSV."
  >>>
  output {
    File stxtyper_report = "~{samplename}_stxtyper.tsv"
    File stxtyper_log = "~{samplename}_stxtyper.log"
    String stxtyper_docker = docker
    String stxtyper_version = read_string("VERSION.txt")
    Int stxtyper_num_hits = read_int("stxtyper_num_hits.txt")
    # outputs parsed from stxtyper output TSV
    String stxtyper_complete_operons = read_string("stxtyper_complete_operons.txt")
    String stxtyper_partial_hits = read_string("stxtyper_partial_hits.txt")
    String stxtyper_stxA_complete_hits = read_string("stxtyper_stxA_complete_hits.txt")
    String stxtyper_stxB_complete_hits = read_string("stxtyper_stxB_complete_hits.txt")
    String stxtyper_stxA_partial_hits = read_string("stxtyper_stxA_partial_hits.txt")
    String stxtyper_stxB_partial_hits = read_string("stxtyper_stxB_partial_hits.txt")
    String stxtyper_stx_frameshifts_or_internal_stop_hits =  read_string("stxtyper_stx_frameshifts_or_internal_stop_hits.txt")
    String stxtyper_A_partial_B_complete = read_string("stxtyper_A_partial_B_complete.txt")
    String stx_novel_hits = read_string("stx_novel_hits.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}
