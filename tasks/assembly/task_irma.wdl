version 1.0

task irma {
  input {
    File read1
    File? read2
    String seq_method
    String samplename
    Boolean keep_ref_deletions = true # set DEL_TYPE config in irma_config.sh to "DEL" if false, "NNN" if true
    Int minimum_consensus_support = 50 # IRMA default is 1, but matching MIRA standards for ONT = 50 and ILMN = 30 via defaults at theiacov workflow level WDLs: https://cdcgov.github.io/MIRA/articles/sequence-qc.html
    Int minimum_read_length = 75 # matching default for TheiaCoV_Illumina_PE; NOTE: IRMA's default is 125 bp
    Int minimum_average_consensus_allele_quality = 10 # IRMA default is 0, we are matching MIRA standards for both ONT and ILMN: https://cdcgov.github.io/MIRA/articles/sequence-qc.html
    Float minimum_ambiguous_threshold = 0.20
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/irma:1.2.0"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 100
  }
  String read_basename = basename(read1)
  command <<<
    # capture irma vesion
    IRMA | head -n1 | awk -F' ' '{ print "IRMA " $5 }' | tee VERSION

    # set -euo pipefail to avoid silent failure; must happen AFTER running "IRMA" since it throws exit code 1
    set -euo pipefail

    ### IRMA configuration ###
    # CPU config
    num_cpus_actual=$(nproc)
    echo "DEBUG: Number of CPUs available: ${num_cpus_actual}. Setting this in irma_config.sh file..."
    echo "SINGLE_LOCAL_PROC=${num_cpus_actual}" > irma_config.sh
    # set this variable to half the value of num_cpus_actual, as per the IRMA documentation: https://wonder.cdc.gov/amd/flu/irma/configuration.html
    echo "DOUBLE_LOCAL_PROC=$((${num_cpus_actual}/2))" >> irma_config.sh

    # any base with less than the minimum support depth will be called N
    echo "MIN_CONS_SUPPORT=~{minimum_consensus_support}" >> irma_config.sh
    # any base with less than the minimum quality will be called N
    echo "MIN_CONS_QUALITY=~{minimum_average_consensus_allele_quality}" >> irma_config.sh

    # minimum called SNV frequency for mixed base calls in amended consensus
    echo "MIN_AMBIG=~{minimum_ambiguous_threshold}" >> irma_config.sh

    # this is done so that IRMA used PWD as the TMP directory instead of /tmp/root that it tries by default; cromwell doesn't allocate much disk space here (64MB or some small amount)
    echo "DEBUG: creating an optional IRMA configuration file to set TMP directory to $(pwd)"
    echo "TMP=$(pwd)" >> irma_config.sh

    # set how to handle deletions
    if ~{keep_ref_deletions}; then 
      echo 'DEL_TYPE="NNN"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    else # default in WDL and IRMA
      # IRMA docs state: If sites are completely missing during read gathering use the reference seed (REF), delete by ambiguation (NNN), or just remove (DEL)
      echo 'DEL_TYPE="DEL"' >> irma_config.sh
      echo 'ALIGN_PROG="BLAT"' >> irma_config.sh
    fi

    # set the minimum read length; matching default for TheiaCoV_Illumina_PE which is 75 bp
    echo 'MIN_LEN=~{minimum_read_length}' >> irma_config.sh

    echo "DEBUG: Custom irma_config.sh file contents:"
    cat irma_config.sh
    echo "DEBUG: End of custom irma_config.sh file contents"
    ### END IRMA CONFIG ###

    # run IRMA
    # set IRMA module depending on sequencing technology
    if [[ ~{seq_method} == "OXFORD_NANOPORE" ]]; then
      IRMA "FLU-minion" "~{read1}" ~{samplename} --external-config irma_config.sh
    else
      # else, assume Illumina paired-end data as input
      IRMA "FLU" "~{read1}" "~{read2}" ~{samplename} --external-config irma_config.sh
    fi

    # capture some IRMA log & config files; rename to use .tsv suffix instead of .txt
    mv -v ~{samplename}/tables/READ_COUNTS.txt ~{samplename}/tables/READ_COUNTS.tsv
    mv -v ~{samplename}/logs/run_info.txt ~{samplename}/logs/run_info.tsv

    # capture IRMA type
    if compgen -G "~{samplename}/*fasta"; then
      # look at list of files that match the above pattern, grab the first one, and extract the type from the filename. We expect: ~{samplename}/B_HA.fasta
      echo "Type_"$(basename "$(echo "$(find ~{samplename}/*.fasta | head -n1)")" | cut -d_ -f1) > IRMA_TYPE
      # set irma_type bash variable which is used later
      irma_type=$(cat IRMA_TYPE)
      
      # flu segments from largest to smallest
      # Thank you Molly H. for this code block!
      # declare associative arrays for segment numbers
      declare -A FluA=(["PB2"]="1" ["PB1"]="2" ["PA"]="3" ["HA"]="4" ["NP"]="5" ["NA"]="6" ["MP"]="7" ["NS"]="8" )
      declare -A FluB=(["PB1"]="1" ["PB2"]="2" ["PA"]="3" ["HA"]="4" ["NP"]="5" ["NA"]="6" ["MP"]="7" ["NS"]="8" )
      # create new array SEGMENT_DICT based on Flu Type (either A or B) to use in the loop below
      if [[ "${irma_type}" == "Type_A" ]]; then
        echo "DEBUG: IRMA type is A. Using FluA array for segment order...."
        declare -A SEGMENT_DICT
        for key in "${!FluA[@]}"; do
          SEGMENT_DICT["$key"]="${FluA["$key"]}"
        done
      elif [[ "${irma_type}" == "Type_B" ]]; then
        echo "DEBUG: IRMA type is B. Using FluB array for segment order...."
        declare -A SEGMENT_DICT
        for key in "${!FluB[@]}"; do
          SEGMENT_DICT["$key"]="${FluB["$key"]}"
        done
      fi

      # for use in below code block as well as at the end for storing padded FASTA files
      # we are keeping this dir OUTSIDE of the IRMA output dir "~{samplename}" so that we can differentiate these files (manually created by task) from the IRMA output files
      mkdir padded_assemblies/

      echo "DEBUG: creating IRMA FASTA file containing all segments in order (largest to smallest)...."
      ### concatenate files in the order of the FluA or FluB segments array ###
      # for each segment number in the SEGMENT_DICT array (sorted by order in FluA/FluB arrays), find the file that contains the segment number in the filename
      for SEGMENT_NUM in $(echo "${SEGMENT_DICT[@]}" | tr ' ' '\n'| sort); do
        segment_file=$(find "~{samplename}/amended_consensus" -name "*_${SEGMENT_NUM}.fa")
        echo "DEBUG: segment_file is set to: $segment_file"
        # if the segment file exists, rename it and added to the final multi FASTA file "~{samplename}.irma.consensus.fasta"
        if [ -n "$segment_file" ]; then
          # craziness so we can use the segment number to get the segment string
          for SEGMENT_STR in "${!SEGMENT_DICT[@]}"; do
            if [[ ${SEGMENT_DICT[$SEGMENT_STR]} == "$SEGMENT_NUM" ]]; then
              echo "DEBUG: File containing ${SEGMENT_STR} found. Now adjusting FASTA header & renaming FASTA file to include segment abbreviation ${SEGMENT_STR}...."
              sed -i "1s|[0-9]$|${SEGMENT_STR}|" "${segment_file}"

              # final filename should be: ~{samplename}/amended_consensus/~{samplename}_HA.fasta
              mv -v "${segment_file}" "~{samplename}/amended_consensus/~{samplename}_${SEGMENT_STR}.fasta"
              # reassign segment_file bash variable to the new filename
              segment_file="~{samplename}/amended_consensus/~{samplename}_${SEGMENT_STR}.fasta"
            fi
          done

          echo "DEBUG: Adding ${segment_file} to consensus FASTA"
          cat "${segment_file}" >> ~{samplename}/amended_consensus/~{samplename}.irma.consensus.fasta
        else
          echo "WARNING: No file containing segment number ${SEGMENT_NUM} found for ~{samplename}"
        fi
      done
      
      echo "DEBUG: creating a smushed and hacky concatenated copy of the IRMA FASTA file..."
      # remove all newlines from the FASTA file to create a single line FASTA file and remove all headers and create a new headerline
      grep -v "^>" ~{samplename}/amended_consensus/~{samplename}.irma.consensus.fasta | tr -d '\n' | sed '1i >~{samplename}_irma_concatenated' > ~{samplename}/amended_consensus/~{samplename}.irma.consensus.concatenated.fasta
      # really hacky way to add an EOF newline but i couldn't be bothered to figure out a better way atm  
      echo "" >> ~{samplename}/amended_consensus/~{samplename}.irma.consensus.concatenated.fasta

      echo "DEBUG: creating copy of consensus FASTA with periods replaced by Ns...."
      # use sed to create copy of FASTA file where periods are replaced by Ns, except in the FASTA header lines that begin with '>'
      sed '/^>/! s/\./N/g' ~{samplename}/amended_consensus/~{samplename}.irma.consensus.fasta > padded_assemblies/~{samplename}.irma.consensus.pad.fasta
    else
      echo "No IRMA assembly generated for flu type prediction" | tee IRMA_TYPE
      echo "Exiting IRMA task early since no IRMA assembly was generated."
      exit 0
    fi

    # rename IRMA outputs to include samplename. Example: "B_HA.fasta" -> "sample0001_HA.fasta"
    echo "DEBUG: Renaming IRMA output VCFs, FASTAs, and BAMs to include samplename...."
    for irma_out in ~{samplename}/*{.vcf,.fasta,.bam}; do
      new_name="~{samplename}_"$(basename "${irma_out}" | cut -d "_" -f2- )
      mv -v "${irma_out}" "${new_name}"
    done
    
    # initialize subtype variable as blank so we can check if it's empty later (for Flu B samples)
    subtype=""

    # capture type A subtype
    if compgen -G "~{samplename}_HA*.fasta"; then # check if HA segment exists
      # NOTE: this if block does not get triggered for Flu B samples, because they do not include a subtype in FASTA filename for HA or NA segment
      if [[ "$(ls ~{samplename}_HA*.fasta)" == *"HA_H"* ]]; then # if so, grab H-type if one is identified in assembly header
        subtype="$(basename ~{samplename}_HA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab H-type from last value in under-score-delimited filename
        # rename HA FASTA file to not include subtype in filename
        echo "DEBUG: renaming HA FASTA file to not include subtype in filename...."
        mv -v ~{samplename}_HA*.fasta ~{samplename}_HA.fasta
      fi
      echo "DEBUG: Running sed to change HA segment FASTA header now..."
      # format HA segment to target output name and rename header to include the samplename. Example FASTA header change: ">B_HA" -> ">sample0001_B_HA"
      sed -i "1s/>/>~{samplename}_/" "~{samplename}_HA.fasta"
    fi

    # if there is a file that matches the pattern AND the file contains an N-type in the header, grab the N-type
    # NOTE: this does not get triggered for Flu B samples, because they do not include a subtype in FASTA filename for HA or NA segment
    if compgen -G "~{samplename}_NA*.fasta" && [[ "$(ls ~{samplename}_NA*.fasta)" == *"NA_N"* ]]; then # check if NA segment exists with an N-type identified in header
      subtype+="$(basename ~{samplename}_NA*.fasta | awk -F _ '{print $NF}' | cut -d. -f1)" # grab N-type from last value in under-score-delimited filename
      echo "DEBUG: subtype is set to: ${subtype}"
      # rename NA FASTA file to not include subtype in filename
      echo "DEBUG: renaming NA FASTA file to not include subtype in filename...."
      mv -v ~{samplename}_NA*.fasta ~{samplename}_NA.fasta

      echo "DEBUG: Running sed to change NA FASTA header now..."
      # format NA segment to target output name and rename header to include the samplename
      sed -i "1s/>/>~{samplename}_/" "~{samplename}_NA.fasta"
    fi

    # if bash variable "subtype" is not empty, write it to a file; 
    # otherwise, write a message indicating no subtype was predicted
    if [ -n "${subtype}" ]; then 
      echo "${subtype}" | tee IRMA_SUBTYPE
    else
      echo "No subtype predicted by IRMA" | tee IRMA_SUBTYPE
    fi

    # if "subtype" is "Type_B" then write a note indicating that IRMA does not differentiate between Victoria and Yamagata lineages
    if [[ "${irma_type}" == "Type_B" ]]; then
      echo "IRMA does not differentiate Victoria and Yamagata Flu B lineages. See abricate_flu_subtype output column" | tee IRMA_SUBTYPE_NOTES
    else
      # create empty file
      touch IRMA_SUBTYPE_NOTES
    fi

    if ls "~{samplename}"_HA?*.bam 1> /dev/null 2>&1; then
      echo "DEBUG: Renaming HA BAM files...."
      for file in "~{samplename}"_HA?*.bam; do
        mv -v "$file" "${file%_HA*.bam}_HA.bam"
      done
    fi
    if ls "~{samplename}"_NA?*.bam 1> /dev/null 2>&1; then
      echo "DEBUG: Renaming NA BAM files...."
      for file in "~{samplename}"_NA?*.bam; do
        mv -v "$file" "${file%_NA*.bam}_NA.bam"
      done
    fi
    
    echo "DEBUG: Creating padded FASTA files for each individual segment FASTA, the multi-FASTA, and the concatenated all-segment FASTA...."
    # this loop looks in the PWD for files ending in .fasta, and creates a copy of the file with periods replaced by Ns and dashes are deleted (since they represent gaps)
    for FASTA in ~{samplename}/amended_consensus/*.fasta; do
      # if the renamed file ends in .fasta; then create copy of the file with periods replaced by Ns and dashes removed in all lines except FASTA headers that begin with '>'
      echo "DEBUG: creating padded FASTA file for ${FASTA}"
      BASENAME=$(basename "${FASTA}")
      sed '/^>/! s/\./N/g' "${FASTA}" > "padded_assemblies/${BASENAME%.fasta}.temp.fasta"
      sed '/^>/! s/-//g' "padded_assemblies/${BASENAME%.fasta}.temp.fasta" > "padded_assemblies/${BASENAME%.fasta}.pad.fasta"
      # clean up temporary FASTA files
      rm "padded_assemblies/${BASENAME%.fasta}.temp.fasta"
    done
  >>>
  output {
    # all of these FASTAs are derived from the amended_consensus/*.fa files produced by IRMA
    File? irma_assembly_fasta = "~{samplename}/amended_consensus/~{samplename}.irma.consensus.fasta"
    File? irma_assembly_fasta_concatenated = "~{samplename}/amended_consensus/~{samplename}.irma.consensus.concatenated.fasta"
    File? seg_ha_assembly = "~{samplename}/amended_consensus/~{samplename}_HA.fasta"
    File? seg_na_assembly = "~{samplename}/amended_consensus/~{samplename}_NA.fasta"
    File? seg_pa_assembly = "~{samplename}/amended_consensus/~{samplename}_PA.fasta"
    File? seg_pb1_assembly = "~{samplename}/amended_consensus/~{samplename}_PB1.fasta"
    File? seg_pb2_assembly = "~{samplename}/amended_consensus/~{samplename}_PB2.fasta"
    File? seg_mp_assembly = "~{samplename}/amended_consensus/~{samplename}_MP.fasta"
    File? seg_np_assembly = "~{samplename}/amended_consensus/~{samplename}_NP.fasta"
    File? seg_ns_assembly = "~{samplename}/amended_consensus/~{samplename}_NS.fasta"
    
    # adding these "padded" assemblies as outputs to be passed to VADR and MAFFT (antiviral substitutions tasks)
    # we may remove these outputs in the future if IRMA code is updated to not output periods in the consensus sequences
    File? irma_assembly_fasta_padded = "padded_assemblies/~{samplename}.irma.consensus.pad.fasta"
    File? irma_assembly_fasta_concatenated_padded = "padded_assemblies/~{samplename}.irma.consensus.concatenated.pad.fasta"
    File? seg_ha_assembly_padded = "padded_assemblies/~{samplename}_HA.pad.fasta"
    File? seg_na_assembly_padded = "padded_assemblies/~{samplename}_NA.pad.fasta"
    File? seg_pa_assembly_padded = "padded_assemblies/~{samplename}_PA.pad.fasta"
    File? seg_pb1_assembly_padded = "padded_assemblies/~{samplename}_PB1.pad.fasta"
    File? seg_pb2_assembly_padded = "padded_assemblies/~{samplename}_PB2.pad.fasta"
    File? seg_mp_assembly_padded = "padded_assemblies/~{samplename}_MP.pad.fasta"
    File? seg_np_assembly_padded = "padded_assemblies/~{samplename}_NP.pad.fasta"
    File? seg_ns_assembly_padded = "padded_assemblies/~{samplename}_NS.pad.fasta"

    String irma_type = read_string("IRMA_TYPE")
    String irma_subtype = read_string("IRMA_SUBTYPE")
    String irma_subtype_notes = read_string("IRMA_SUBTYPE_NOTES")
    Array[File] irma_plurality_consensus_assemblies = glob("~{samplename}*.fasta")
    Array[File] irma_vcfs = glob("~{samplename}*.vcf")
    Array[File] irma_bams = glob("~{samplename}*.bam")
    String irma_docker = docker
    String irma_version = read_string("VERSION")
    # tracking this for outputting min depth threshold used for calling consensus nucleotides in IRMA
    Int irma_minimum_consensus_support = minimum_consensus_support
    File irma_read_counts_tsv = "~{samplename}/tables/READ_COUNTS.tsv"
    File irma_run_info_tsv = "~{samplename}/logs/run_info.tsv"
    File irma_nr_read_counts = "~{samplename}/logs/NR_COUNTS_log.txt"
    # for now just adding bams for these segments for mean coverage calculation
    File? seg_ha_bam = "~{samplename}_HA.bam"
    File? seg_na_bam = "~{samplename}_NA.bam"
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