version 1.0

task stxtyper {
  input {
    File assembly
    String samplename
    # this threshold is for display/outputing partial hits to Terra data table; it does NOT change the coverage threshold used by StxTyper to call hits.
    Float coverage_threshold_to_output_partial_hits = 60.0
    String docker = "kapsakcj/stxtyper:78754d7"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 2
  }
  command <<<
    # capture date
    date | tee DATE
    
    # fail task if any commands below fail (AGH!)
    set -eo pipefail

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
      touch stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stxA_complete_hits.txt stxtyper_stxB_complete_hits.txt stxtyper_stxA_partial_hits.txt stxtyper_stxB_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stxtyper_A_partial_B_complete.txt stx_novel_hits.txt stxtyper_A_B_complete_different_contigs.txt
      # put "none" into all of them so task does not fail
      echo "None" | tee stxtyper_complete_operons.txt stxtyper_partial_hits.txt stxtyper_stxA_complete_hits.txt stxtyper_stxB_complete_hits.txt stxtyper_stxA_partial_hits.txt stxtyper_stxB_partial_hits.txt stxtyper_stx_frameshifts_or_internal_stop_hits.txt stxtyper_A_partial_B_complete.txt stx_novel_hits.txt stxtyper_A_B_complete_different_contigs.txt
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

      ### stxA complete hits
      echo "DEBUG: Parsing stxA complete hits..."
      # explanation: if "A_reference_suptype" col contains "stxA" AND %ID >= 98.00 AND %COV >= 100.00, print the "A_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$10 ~ "stxA.*" && $11 >= 98.00 && $12 >= 100.00 {print $10}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxA_complete_hits.txt
      # if no stxA complete hits found, write "None" to file for output string
      if [ "$(grep --silent 'stxA' stxtyper_stxA_complete_hits.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stxA_complete_hits.txt
      fi

      ### stxB complete hits
      echo "DEBUG: Parsing stxB complete hits..."
      # explanation: if "B_reference_suptype" col contains stxB AND %ID >= 98.00 AND %COV >= 100.00, print the "B_reference_suptype" in a sorted/unique list
      awk -F'\t' -v OFS=, '$14 ~ "stxB.*" && $15 >= 98.00 && $16 >= 100.00 {print $14}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxB_complete_hits.txt
      # if no stxB complete hits found, write "None" to file for output string
      if [ "$(grep --silent 'stx' stxtyper_stxB_complete_hits.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stxB_complete_hits.txt
      fi

      ### partial hits (to any gene in stx operon)
      echo "DEBUG: Parsing stxtyper partial hits..."
      # explanation: if "operon" column contains "PARTIAL" (either PARTIAL or PARTIAL_CONTIG_END possible); print either "stx1" or "stx2" or "stx1,stx2"
      awk -F'\t' -v OFS=, '$4 ~ "PARTIAL.*" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_partial_hits.txt
      # if no stx partial hits found, write "None" to file for output string
      if [ "$(grep --silent 'stx' stxtyper_partial_hits.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_partial_hits.txt
      fi

      ### stxA partial hits
      echo "DEBUG: Parsing stxA partial hits..."
      # explanation: if "A_reference_suptype" col contains "stxA" AND %COV > coverage_threshold_to_output_partial_hits AND %COV < 100.00, print the "A_reference_suptype" in a sorted/unique list
      awk -F'\t' -v cov_threshold=~{coverage_threshold_to_output_partial_hits} -v OFS=, '$10 ~ "stxA.*" && $12 >= cov_threshold && $12 < 100.00 {print $10}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stxA_partial_hits.txt
      # if no stxA partial hits found, write "None" to file for output string
      if [ "$(grep --silent 'stx' stxtyper_stxA_partial_hits.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stxA_partial_hits.txt
      fi

      #### COMMENTING OUT AT USER'S REQUEST, KEEPING CODE IN CASE WE WANT TO BRING BACK LATER ###
      #### stxB partial hits
      #echo "DEBUG: Parsing stxB partial hits..."
      ## explanation: if "B_reference_suptype" col contains stxB AND %COV > coverage_threshold_to_output_partial_hits AND %COV < 100.00, print the "B_reference_suptype" in a sorted/unique list
      #awk -F'\t' -v cov_threshold=~{coverage_threshold_to_output_partial_hits} -v OFS=, '$14 ~ "stxB.*" && $16 >= cov_threshold && $16 < 100.00 {print $14}' ~{samplename}_stxtyper.tsv | sort | #uniq | paste -sd, - | tee stxtyper_stxB_partial_hits.txt
      ## if no stxB partial hits found, write "None" to file for output string. If grep fails to find "stx", check exit code and if greater than 0, write "None" to file
      #if [ "$(grep --silent 'stx' stxtyper_stxB_partial_hits.txt; echo $?)" -gt 0 ]; then
      #  echo "None" > stxtyper_stxB_partial_hits.txt
      #fi

      ### frameshifts or internal stop codons in stx genes
      echo "DEBUG: Parsing stx frameshifts or internal stop codons..."
      # explanation: if operon column contains "FRAME_SHIFT" or "INTERNAL_STOP", print the "operon" in a sorted/unique list
      awk -F'\t' -v OFS=, '$4 == "FRAMESHIFT" || $4 == "INTERNAL_STOP" {print $3}' ~{samplename}_stxtyper.tsv | sort | uniq | paste -sd, - | tee stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      # if no frameshifts or internal stop codons found, write "None" to file for output string
      if [ "$(grep --silent -E 'FRAMESHIFT|INTERNAL_STOP' ~{samplename}_stxtyper.tsv; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_stx_frameshifts_or_internal_stop_hits.txt
      fi

      ### A subunit is partial, but B subunits is complete; can be on same contig OR different contigs
      echo "DEBUG: Parsing A partial, B complete hits of the same subtype..."
      # first create line-delimited list of partial and complete lists of stxA and stxB hits (instead of comma deliminted) for easier parsing
      tr ',' '\n' < stxtyper_stxA_partial_hits.txt > stxtyper_stxA_partial_hits_lines.txt
      tr ',' '\n' < stxtyper_stxB_complete_hits.txt > stxtyper_stxB_complete_hits_lines.txt

      # then check for A partial, B complete hits of the same subtype
      # sed to strip out stxA or stxB from subtype names; then look for matching B complete hits
      echo "DEBUG: searching for subtype matches between stxA partial and stxB complete hits..."
      sed 's/stx[A,B]//g' stxtyper_stxA_partial_hits_lines.txt | while read -r SUBTYPE; do 
        if [ "$(grep --silent "$SUBTYPE" stxtyper_stxB_complete_hits_lines.txt; echo $?)" -eq 0 ]; then
          grep "$SUBTYPE" stxtyper_stxB_complete_hits_lines.txt | paste -sd, - | tee stxtyper_A_partial_B_complete_list.txt; 
        else
          echo "DEBUG: No matching B complete hits found for $SUBTYPE"
        fi
      done
      # create empty file in case there are no subtype matches between stxA partial and stxB complete hits ( so that sed cmd 4 lines down doesn't fail)
      touch stxtyper_A_partial_B_complete_list.txt
      # strip out captial A or B so we are left with list of subtypes, example: stx2c (and NOT stxB2c which describes the B subunit)
      echo "DEBUG: stripping out A or B from gene names to be left with list of subytypes..."
      sed 's/[A,B]//g' stxtyper_A_partial_B_complete_list.txt | paste -sd, - | tee stxtyper_A_partial_B_complete.txt
      # if list is empty, write "None" to file
      if [ "$(grep --silent 'stx' stxtyper_A_partial_B_complete.txt; echo $?)" -gt 0 ]; then
        echo "None" > stxtyper_A_partial_B_complete.txt
      fi

      ### A and B subunits are complete, but on different contigs
      echo "DEBUG: Parsing A and B complete hits on different contigs..."
      # first create line-delimited list of A complete hits of stxA (stxB list made earlier) for easier parsing
      tr ',' '\n' < stxtyper_stxA_complete_hits.txt > stxtyper_stxA_complete_hits_lines.txt
      # create list of stxA with >98.0%ID and 100%COV and include contig name
      echo 'DEBUG: creating list of stxA with >98.0%ID and 100%COV; including contig name...'
      awk -F'\t' -v OFS=, '$10 ~ "stxA.*" && $11 >= 98.00 && $12 >= 100.00 {print $10,$2}' ~{samplename}_stxtyper.tsv | tee stxtyper_stxA_complete_hits_with_contig.txt
      # create list of stxB with >98.0%ID and 100%COV and include contig name
      echo 'DEBUG: creating list of stxB with >98.0%ID and 100%COV; including contig name...'
      awk -F'\t' -v OFS=, '$14 ~ "stxB.*" && $15 >= 98.00 && $16 >= 100.00 {print $14,$2}' ~{samplename}_stxtyper.tsv | tee stxtyper_stxB_complete_hits_with_contig.txt
      # create empty file in case there are no subtype matches between stxB complete hits
      touch stxB_subtype-matched_contig_name.txt
      
     # for every stxA complete hit, search for stxB complete hits of the SAME subtype, then compare contig names
       # STEPS:
       # 1 loop through stxtyper_stxA_complete_hits_with_contig.txt
       #   2 grep for the SUBTYPE in stxtyper_stxB_complete_hits_with_contig.txt
       #   3 if grep is successful (meaning there is a matching stxB hit with same subtype), 
       #      3a compare contig names. if they match, do nothing; if they do NOT match, write the hits (including the contig name) to output txt file (meaning there is complete stxA and stxB hits on diff contigs!)
       #   4 if grep FAILS (meaning there is NOT a stxB hit with the same subtype), then do nothing
       #   5 Last check - do we have complete stxA and stxB hits on diff contigs?
       #     grep for 'stx' in final output txt file, if successful do nothing; if failure (because it's blank), write "None" to output txt file for final output

      # compare 2 lists to see if any stxA and stxB (COMPLETE) of the same subtype are on different contigs
      # 1 loop through stxtyper_stxA_complete_hits_with_contig.txt. Example $HIT_AND_CONTIG: "stxA2c,contig00178"
      cat stxtyper_stxA_complete_hits_with_contig.txt | while read -r HIT_AND_CONTIG; do 
        HIT=$(echo "$HIT_AND_CONTIG" | cut -d ',' -f1) # stxA2c
        echo "DEBUG: HIT is set to: $HIT"
        SUBTYPE=$(echo "$HIT_AND_CONTIG" | sed 's/stx[A,B]//g; s/,.*//g') # 2c
        echo "DEBUG: SUBTYPE is set to: $SUBTYPE"
        CONTIG=$(echo "$HIT_AND_CONTIG" | cut -d ',' -f2) # contig00178
        echo "DEBUG: CONTIG is set to: $CONTIG"
        # 2 grep for the SUBTYPE in stxtyper_stxB_complete_hits_with_contig.txt
        if [ "$(grep --silent "$SUBTYPE" stxtyper_stxB_complete_hits_with_contig.txt; echo $?)" -eq 0 ]; then 
          # 3 if grep is successful (meaning there is a matching stxB hit with same subtype),
          # 3a compare contig names. if they match, do nothing; if they do NOT match, write the hits (including the contig name) to output txt file (meaning there is complete stxA and stxB hits on diff contigs!)
          if [ "$(grep --silent "${CONTIG}" stxtyper_stxB_complete_hits_with_contig.txt; echo $?)" -gt 0 ]; then
            B_HIT=$(grep "$SUBTYPE" stxtyper_stxB_complete_hits_with_contig.txt | cut -d ',' -f1)
            echo "DEBUG: B_HIT is set to: $B_HIT"
            B_CONTIG=$(grep "$SUBTYPE" stxtyper_stxB_complete_hits_with_contig.txt | cut -d ',' -f2)
            echo "DEBUG: B_CONTIG is set to: $B_CONTIG"
            # if contig names do NOT match, write the hits (including the contig name) to output txt file
            echo "${HIT},${CONTIG},${B_HIT},${B_CONTIG};" >> stxtyper_A_B_complete_different_contigs.txt
          else
            B_HIT=$(grep "$SUBTYPE" stxtyper_stxB_complete_hits_with_contig.txt | cut -d ',' -f1)
            echo "DEBUG: B_HIT is set to: $B_HIT"
            echo "DEBUG: Contig names match for $HIT and ${B_HIT}, not adding to list of complete hits on different contigs."
          fi
        fi
      done
      # now that this mess of a code block has finished running, check to see if it's blank (meaning no complete stxA and stxB of the same subtype found on diff contigs), if it's blank, write "None" to file
      # touch the file in case it doesn't exist yet (also meaning no complete stxA and stxB of the same subtype found on diff contigs)
      touch stxtyper_A_B_complete_different_contigs.txt
      if [ "$(grep --silent 'stx' stxtyper_A_B_complete_different_contigs.txt; echo $?)" -gt 0 ]; then
        echo "DEBUG: No complete stxA and stxB hits of the same subtype found on different contigs. Writing 'None' to output string"
        echo "None" > stxtyper_A_B_complete_different_contigs.txt
      fi
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
    # commenting out at user's request. keeping code in case we want to bring back later
    #String stxtyper_stxB_partial_hits = read_string("stxtyper_stxB_partial_hits.txt")
    String stxtyper_stx_frameshifts_or_internal_stop_hits =  read_string("stxtyper_stx_frameshifts_or_internal_stop_hits.txt")
    String stxtyper_A_partial_B_complete = read_string("stxtyper_A_partial_B_complete.txt")
    String stxtyper_A_B_complete_different_contigs = read_string("stxtyper_A_B_complete_different_contigs.txt")
    String stxtyper_novel_hits = read_string("stx_novel_hits.txt")
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
