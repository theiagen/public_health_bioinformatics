version 1.0

task vadr_flu_segments {
  input {
    File genome_fasta
    File vadr_outputs_tgz
    Int cpu = 1
    Int memory = 2
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/ubuntu/ubuntu:jammy-20230816"
  }
  String out_base = basename(vadr_outputs_tgz, ".vadr.tar.gz")
  command <<<
    set -euo pipefail

    # extract the tarball
    mkdir -p ~{out_base}
    tar -C ~{out_base} -xzf ~{vadr_outputs_tgz}

    # skip header lines and convert classification file to tsv format (21 fields)
    tail -n +4 ~{out_base}/~{out_base}.vadr.sqc | sed -E 's/  +/\t/g' > ~{out_base}/~{out_base}.vadr.sqc.tsv

    # does the assembly contain influenza segments?
    grp1_model=$(awk '/flu[AB]-seg/' ~{out_base}/~{out_base}.vadr.sqc.tsv)

    if [[ -n "$grp1_model" ]]; then
      # split the assembly fasta into segments based on headers
      mkdir -p tmp_fastas/
      csplit -s -z -f tmp_fastas/seq_ -b "%d" ~{genome_fasta} '/^>/' '{*}'

      # create final directory for flu segments
      mkdir -p flu_segments/

      # define segment mappings using associative arrays. derived from: https://github.com/ncbi/vadr/wiki/Influenza-annotation#influenza-vadr-model-library-1
      declare -A fluA_map=( ["seg1"]="PB2" ["seg2"]="PB1" ["seg3"]="PA" ["seg4"]="HA" ["seg5"]="NP" ["seg6"]="NA" ["seg7"]="MP" ["seg8"]="NS" )
      declare -A fluB_map=( ["seg1"]="PB1" ["seg2"]="PB2" ["seg3"]="PA" ["seg4"]="HA" ["seg5"]="NP" ["seg6"]="NA" ["seg7"]="MP" ["seg8"]="NS" )

      # loop through VADR classification summary and extract the available flu type/gene segments
      while IFS=$'\t' read -r -a line; do
        SEQ_NAME="${line[1]}"
        MODEL="${line[6]}"
        FLU_TYPE=$(echo "${MODEL}" | cut -d'-' -f1)
        NUM_SEGMENT=$(echo "${MODEL}" | cut -d'-' -f2)

        # select map based on flu type
        if [[ "${FLU_TYPE}" == "fluA" ]]; then
          GENE_SEGMENT=${fluA_map["${NUM_SEGMENT}"]}
        else
          GENE_SEGMENT=${fluB_map["${NUM_SEGMENT}"]}
        fi

        echo "SEQ_NAME: ${SEQ_NAME}, MODEL: ${MODEL}, FLU_TYPE: ${FLU_TYPE}, NUM_SEGMENT: ${NUM_SEGMENT}, GENE_SEGMENT: ${GENE_SEGMENT}"

        # rename the segment fasta file if the header matches the sequence name (pulled from the VADR classification summary)
        for fasta in tmp_fastas/seq_*; do
          # check first line if it starts with '>', get only the first field without the '>'
          HEADER_NAME=$(awk 'NR==1 && /^>/ {sub(/^>/, ""); print $1}' "${fasta}")

          if [[ "${HEADER_NAME}" == "${SEQ_NAME}" ]]; then
            output_segment_fasta="flu_segments/~{out_base}_${GENE_SEGMENT}.fasta"
            cp "$fasta" "$output_segment_fasta"
          fi
        done
      done < ~{out_base}/~{out_base}.vadr.sqc.tsv

      # defining the concatenated fasta header
      CONCATENATED_HEADER=">~{out_base}_concatenated_flu_segments"

      # concatenate all segments into a single fasta file in the correct order (seg1 -> seg8)
      # need a separate loop because the VADR classification summary does not always report segments in the correct order
      for i in {1..8}; do
        NUM_SEGMENT="seg${i}"

        # select map based on the last flu type identified by VADR. note that flu A and B cannot recombine, so the last reported flu type should be consistent across all segments
        if [[ "${FLU_TYPE}" == "fluA" ]]; then
          GENE_SEGMENT=${fluA_map["${NUM_SEGMENT}"]}
        else
          GENE_SEGMENT=${fluB_map["${NUM_SEGMENT}"]}
        fi

        SEGMENT_FASTA_FILE="flu_segments/~{out_base}_${GENE_SEGMENT}.fasta"
        if [[ -s "$SEGMENT_FASTA_FILE" ]]; then
          # append gene segment to the concatenated fasta header
          CONCATENATED_HEADER="${CONCATENATED_HEADER}-${GENE_SEGMENT}"

          # append sequence and remove header/newlines
          grep -v ">" "$SEGMENT_FASTA_FILE" | tr -d '\n' >> "~{out_base}_concatenated.fasta"
        else
          echo "Warning: Missing segment fasta file $SEGMENT_FASTA_FILE. Does not exist or is empty. Skipping concatenation for this segment."
        fi
      done

      # Insert the header at the beginning of the file
      sed -i "1i${CONCATENATED_HEADER}" "~{out_base}_concatenated.fasta"

      # add final newline to the concatenated fasta
      echo "" >> "~{out_base}_concatenated.fasta"
    fi
  >>>
  output {
    File? assembly_fasta_concatenated = "~{out_base}_concatenated.fasta"
    File? seg_ha_assembly = "flu_segments/~{out_base}_HA.fasta"
    File? seg_na_assembly = "flu_segments/~{out_base}_NA.fasta"
    File? seg_pa_assembly = "flu_segments/~{out_base}_PA.fasta"
    File? seg_pb1_assembly = "flu_segments/~{out_base}_PB1.fasta"
    File? seg_pb2_assembly = "flu_segments/~{out_base}_PB2.fasta"
    File? seg_mp_assembly = "flu_segments/~{out_base}_MP.fasta"
    File? seg_np_assembly = "flu_segments/~{out_base}_NP.fasta"
    File? seg_ns_assembly = "flu_segments/~{out_base}_NS.fasta"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
  }
}