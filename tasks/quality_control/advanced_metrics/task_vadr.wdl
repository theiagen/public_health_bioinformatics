version 1.0

task vadr {
  meta {
    description: "Runs NCBI's Viral Annotation DefineR for annotation and QC. See https://github.com/ncbi/vadr/wiki/Coronavirus-annotation"
  }
  input {
    File genome_fasta
    String vadr_opts = "--noseqnamemax --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --out_allfasta"
    Int assembly_length_unambiguous
    Int skip_length = 10000
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.6.3-hav-flu2"
    Int min_length = 50
    Int max_length = 30000
    Int cpu = 4
    Int memory = 16
    Int disk_size = 100
  }
  String out_base = basename(basename(basename(genome_fasta, ".fasta"), ".fa"), ".fna")
  command <<<
    set -euo pipefail

    # set default values for flu type and subtype
    echo "N/A" > FLU_SUBTYPE
    echo "N/A" > FLU_TYPE

    if [ ~{assembly_length_unambiguous} -gt ~{skip_length} ]; then

      # remove terminal ambiguous nucleotides
      /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
        "~{genome_fasta}" \
        --minlen ~{min_length} \
        --maxlen ~{max_length} \
        > "~{out_base}_trimmed.fasta"

      # run VADR
      # --split and --cpu must be used in conjuction
      v-annotate.pl \
        --split --cpu ~{cpu} \
        ~{vadr_opts} \
        "~{out_base}_trimmed.fasta" \
        "~{out_base}"

      # package everything for output
      tar -C "~{out_base}" -czvf "~{out_base}.vadr.tar.gz" .

      # skip header lines and convert classification file to tsv format (21 fields)
      tail -n +4 ~{out_base}/~{out_base}.vadr.sqc | sed -E 's/  +/\t/g' > ~{out_base}/~{out_base}.vadr.sqc.tsv

      # does the assembly contain influenza segments?
      grp1_model=$(awk '/PASS/ && /flu[AB]-seg/' ~{out_base}/~{out_base}.vadr.sqc.tsv)

      if [[ -n "$grp1_model" ]]; then
        # split the assembly fasta into segments based on headers
        mkdir -p tmp_fastas/
        csplit -s -z -f tmp_fastas/seq_ -b "%d" ~{genome_fasta} '/^>/' '{*}'

        # create final directory for flu segments
        mkdir -p flu_segments/

        # define segment mappings using associative arrays
        declare -A fluA_map=( ["seg1"]="PB2" ["seg2"]="PB1" ["seg3"]="PA" ["seg4"]="HA" ["seg5"]="NP" ["seg6"]="NA" ["seg7"]="MP" ["seg8"]="NS" )
        declare -A fluB_map=( ["seg1"]="PB1" ["seg2"]="PB2" ["seg3"]="PA" ["seg4"]="HA" ["seg5"]="NP" ["seg6"]="NA" ["seg7"]="MP" ["seg8"]="NS" )

        HA_SUBTYPE=""
        NA_SUBTYPE=""
        # extract the flu type, subtype and gene segment from the classification summary
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

          # grab HA/NA subtype if available
          if [[ "${line[7]}" != "-" && "${MODEL}" == "fluA-seg4" ]]; then
            HA_SUBTYPE="${line[7]}"
          fi
          if [[ "${line[7]}" != "-" && "${MODEL}" == "fluA-seg6" ]]; then
            NA_SUBTYPE="${line[7]}"
          fi

          # create concatentated fasta file header
          echo ">vadr_concatenated_flu_segments" > "~{out_base}_concatenated.fasta"

          # find matching fasta file and copy it
          for fasta in tmp_fastas/seq_*; do
            # check if the sequence name is in the header line of the fasta file
            if head -n 1 "${fasta}" | grep "${SEQ_NAME}"; then
              output_segment_fasta="flu_segments/~{out_base}_${GENE_SEGMENT}.fasta"
              cp "$fasta" "$output_segment_fasta"
              grep -v ">" "$fasta" >> "~{out_base}_concatenated.fasta"
              break
            fi
          done
        done < ~{out_base}/~{out_base}.vadr.sqc.tsv

        FLU_SUBTYPE="${HA_SUBTYPE}${NA_SUBTYPE}"
        echo "${FLU_SUBTYPE}" > FLU_SUBTYPE
        echo "${FLU_TYPE}" > FLU_TYPE

      fi

      # package up FASTA files into zip file for output. Note: this will work whether the --out_allfasta flag is included or not (there are just more when the option is used)
      mkdir -v vadr_fasta_files
      cp -v ~{out_base}/*.fa vadr_fasta_files
      zip ~{out_base}_vadr-fasta-files.zip vadr_fasta_files/*.fa 

      # prep alerts into a tsv file for parsing
      cut -f 5 "~{out_base}/~{out_base}.vadr.alt.list" | tail -n +2 > "~{out_base}.vadr.alerts.tsv"
      cat "~{out_base}.vadr.alerts.tsv" | wc -l > NUM_ALERTS

      # rename sequence classification summary file to end in txt
      mv -v "~{out_base}/~{out_base}.vadr.sqc" "~{out_base}/~{out_base}.vadr.sqc.txt"

    else
      echo "VADR skipped due to poor assembly; assembly length (unambiguous) = ~{assembly_length_unambiguous}" > NUM_ALERTS
    fi

  >>>
  output {
    File? feature_tbl_pass = "~{out_base}/~{out_base}.vadr.pass.tbl"
    File? feature_tbl_fail = "~{out_base}/~{out_base}.vadr.fail.tbl"
    File? classification_summary_file = "~{out_base}/~{out_base}.vadr.sqc.txt"
    String num_alerts = read_string("NUM_ALERTS")
    File? alerts_list = "~{out_base}/~{out_base}.vadr.alt.list"
    File? outputs_tgz = "~{out_base}.vadr.tar.gz"
    File? vadr_fastas_zip_archive = "~{out_base}_vadr-fasta-files.zip"
    String vadr_docker = docker
    String flu_type = read_string("FLU_TYPE")
    String flu_subtype = read_string("FLU_SUBTYPE")
    File? segmented_assemblies_concatenated = "~{out_base}_concatenated.fasta"
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
    memory: memory + " GB"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    cpu: cpu
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}