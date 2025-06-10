version 1.0

task stxtyper {
  input {
    File assembly
    String samplename
    Boolean enable_debugging = false # Additional messages are printed and files in $TMPDIR are not removed after running
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/stxtyper:1.0.42"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 4
  }
  command <<<
    # fail task if any commands below fail since there's lots of bash conditionals below (AGH!)
    set -eo pipefail

    # capture version info
    stxtyper --version | tee VERSION.txt

    # NOTE: by default stxyper uses $TMPDIR or /tmp, so if we run into issues we may need to adjust in the future. Could potentially use PWD as the TMPDIR.
    echo "DEBUG: TMPDIR is set to: $TMPDIR"

    echo "DEBUG: running StxTyper now..."
    # run StxTyper on assembly; may need to add/remove options in the future if they change
    # NOTE: stxtyper can accept gzipped assemblies, so no need to unzip
    stxtyper \
      --nucleotide ~{assembly} \
      --name ~{samplename} \
      --output ~{samplename}_stxtyper.tsv \
      --threads ~{cpu} \
      ~{true='--debug' false='' enable_debugging} \
      --log ~{samplename}_stxtyper.log

    # parse output TSV
    echo "DEBUG: Parsing StxTyper output TSV..."

    # check for output file with only 1 line (meaning no hits found); exit cleanly if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -eq 1 ]; then
      echo "No hits found by StxTyper" > stxtyper_hits.txt
      echo "0" > stxtyper_num_hits.txt
      echo "DEBUG: No hits found in StxTyper output TSV. Exiting task with exit code 0 now."
      
      # create empty output files
      touch stxtyper_all_hits.txt stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt  stx_novel_hits.txt stxtyper_extended_operons.txt stxtyper_ambiguous_hits.txt
      # put "none" into all of them so task does not fail
      echo "None" | tee stxtyper_all_hits.txt stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stx_novel_hits.txt stxtyper_extended_operons.txt stxtyper_ambiguous_hits.txt
      exit 0
    fi
    
    # check for output file with more than 1 line (meaning hits found); count lines & parse output TSV if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -gt 1 ]; then
      echo "Hits found by StxTyper. Counting lines & parsing output TSV now..."
      # count number of lines in output TSV (excluding header)
      wc -l < ~{samplename}_stxtyper.tsv | awk '{print $1-1}' > stxtyper_num_hits.txt
      # remove header line
      sed '1d' ~{samplename}_stxtyper.tsv > ~{samplename}_stxtyper_noheader.tsv

      ##### parse output TSV #####
      ### complete operons
      echo "DEBUG: Parsing complete operons..."
      awk -F'\t' -v OFS=, '$4 == "COMPLETE" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stxtyper_complete_operons.txt
      # if grep for COMPLETE fails, write "None" to file for output string
      if [[ "$(grep --silent 'COMPLETE' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]]; then
        echo "None" > stxtyper_complete_operons.txt
      fi

      ### complete_novel operons
      echo "DEBUG: Parsing complete novel hits..."
      awk -F'\t' -v OFS=, '$4 == "COMPLETE_NOVEL" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stx_novel_hits.txt
      # if grep for COMPLETE_NOVEL fails, write "None" to file for output string
      if [ "$(grep --silent 'COMPLETE_NOVEL' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stx_novel_hits.txt
      fi

      ### partial hits (to any gene in stx operon)
      echo "DEBUG: Parsing stxtyper partial hits..."
      # explanation: if "operon" column contains "PARTIAL" (either PARTIAL or PARTIAL_CONTIG_END possible); print either "stx1" or "stx2" or "stx1,stx2"
      awk -F'\t' -v OFS=, '$4 ~ "PARTIAL.*" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_partial_hits.txt
      # if no stx partial hits found, write "None" to file for output string
      if [ "$(grep --silent 'stx' stxtyper_partial_hits.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_partial_hits.txt
      fi

      ### frameshifts or internal stop codons in stx genes
      echo "DEBUG: Parsing stx frameshifts or internal stop codons..."
      # explanation: if operon column contains "FRAME_SHIFT" or "INTERNAL_STOP", print the "operon" in a sorted/unique list
      awk -F'\t' -v OFS=, '$4 == "FRAMESHIFT" || $4 == "INTERNAL_STOP" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      # if no frameshifts or internal stop codons found, write "None" to file for output string
      if [ "$(grep --silent -E 'FRAMESHIFT|INTERNAL_STOP' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      fi

      ### extended operons
      echo "DEBUG: Parsing extended operons..."
      awk -F'\t' -v OFS=, '$4 == "EXTENDED" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stxtyper_extended_operons.txt
      if [ "$(grep --silent 'EXTENDED' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_extended_operons.txt
      fi

      ### ambiguous hits
      echo "DEBUG: Parsing ambiguous hits..."
      awk -F'\t' -v OFS=, '$4 == "AMBIGUOUS" {print $3}' ~{samplename}_stxtyper.tsv | paste -sd, - | tee stxtyper_ambiguous_hits.txt
      if [ "$(grep --silent 'AMBIGUOUS' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_ambiguous_hits.txt
      fi
      
      echo "DEBUG: generating stx_type_all string output now..."
      # sort and uniq so there are no duplicates; then paste into a single comma-separated line with commas
      # sed is to remove any instances of "None" from the output
      cat stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stx_novel_hits.txt stxtyper_extended_operons.txt stxtyper_ambiguous_hits.txt | sed '/None/d' | sort | uniq | paste -sd, - > stxtyper_all_hits.txt

    fi
    echo "DEBUG: Finished parsing StxTyper output TSV."
  >>>
  output {
    File stxtyper_report = "~{samplename}_stxtyper.tsv"
    File stxtyper_log = "~{samplename}_stxtyper.log"
    String stxtyper_docker = docker
    String stxtyper_version = read_string("VERSION.txt")
    # outputs parsed from stxtyper output TSV
    Int stxtyper_num_hits = read_int("stxtyper_num_hits.txt")
    String stxtyper_all_hits = read_string("stxtyper_all_hits.txt")
    String stxtyper_complete_operon_hits = read_string("stxtyper_complete_operons.txt")
    String stxtyper_partial_hits = read_string("stxtyper_partial_hits.txt")
    String stxtyper_frameshifts_or_internal_stop_hits =  read_string("stxtyper_stx_frameshifts_or_internal_stop_hits.txt")
    String stxtyper_novel_hits = read_string("stx_novel_hits.txt")
    String stxtyper_extended_operons = read_string("stxtyper_extended_operons.txt")
    String stxtyper_ambiguous_hits = read_string("stxtyper_ambiguous_hits.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1 # does not take long (usually <3 min) to run stxtyper on 1 genome, preemptible is fine
    maxRetries: 3
  }
}
