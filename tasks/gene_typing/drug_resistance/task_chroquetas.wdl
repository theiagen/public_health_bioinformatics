version 1.0

task chroquetas {
  input {
    File assembly_fasta
    String species
    String samplename

    Int min_percent_coverage = 75 # set to ChroQueTas default
    Int min_percent_identity = 90  # set to ChroQueTas default
    File? proteome_fasta
    String? translation_code

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/chroquetas:1.0.0"
    Int cpu = 2
    Int disk_size = 16
    Int memory = 8
  }
  command <<<
  # Multiple piped commands return non-zero exit status, so -o is incompatible
  set -eu

  # Monitor status with CHROQUETAS_STATUS file
  echo "ERROR" > CHROQUETAS_STATUS

  # get version (non-zero exit status)
  ChroQueTas.sh --version \
    | sed -E 's/\(Chromosome Query Targets\) version //' \
    | tee CHROQUETAS_VERSION

  # ChroQueTas expects species with underscore delimiters and a capital first letter
  species_prep=$(echo ~{species} | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
  corrected_species=$(echo ${species_prep^})

  # ensure the species is compatible
  grep_check=$(ChroQueTas.sh --list_species | cut -f 1 | tail -n+2 | grep -P "^${corrected_species}$" || true)
  if [ -z ${grep_check} ]; then
    echo "ERROR: Incompatible species" | tee CHROQUETAS_STATUS
    echo "" > AMR_SUMMARY_STRING
    echo "" > ANNOTATED_AMR_SUMMARY_STRING
  else
    # call chroquetas
    ChroQueTas.sh \
      -g ~{assembly_fasta} \
      -s ${corrected_species} \
      --min_cov ~{min_percent_coverage} \
      --min_id ~{min_percent_identity} \
      -t ~{cpu} \
      -o chroquetas_out \
      ~{if defined(proteome_fasta) then "-p ~{proteome_fasta}" else ""} \
      ~{if defined(translation_code) then "-c ~{translation_code}" else ""}
  
    # rename output files to include sample name
    based_name=$(echo $(basename ~{assembly_fasta}) | sed -E 's/(.*)\.[^\.]+$/\1/')
    if [ -f chroquetas_out/${based_name}.ChroQueTaS.AMR_stats.txt ]; then
      mv chroquetas_out/${based_name}.ChroQueTaS.AMR_stats.txt chroquetas_out/~{samplename}.ChroQueTaS.AMR_stats.txt
      mv chroquetas_out/${based_name}.ChroQueTaS.AMR_summary.txt chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt
    
      # extract AMR summary string
      # e.g. <GENE>_<REF_POSITION><AA_CHANGE><QUERY_POSITION>
      tail -n+2 chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt \
        | awk '{ print $1, $4, $3, $5 }' \
        | sed -E 's/([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)/\1_\2\3\4/' \
        | tr '\n' ',' \
        | sed -E 's/,$//' \
        | tee AMR_SUMMARY_STRING
  
      # extract AMR summary annotated w/fungicide resistance
      # e.g. <GENE>_<REF_POSITION><AA_CHANGE><QUERY_POSITION>[<FUNGICIDE_RESISTANCE1>;<FUNGICIDE_RESISTANCEn>]
      tail -n+2 chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt \
        | awk '{ print $1, $4, $3, $5, $6 }' \
        | sed -E 's/([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)/\1_\2\3\4[\5]/' \
        | sed -E 's/,/;/g' \
        | tr '\n' ',' \
        | sed -E 's/,$//' \
        | tee ANNOTATED_AMR_SUMMARY_STRING
  
      echo "PASS" | tee CHROQUETAS_STATUS
    # check if no AMR genes found
    elif [ -f chroquetas_out/${based_name}.ChroQueTaS.AMR_summary.txt ]; then
      # if no AMR genes found, create empty output files
      mv chroquetas_out/${based_name}.ChroQueTaS.AMR_summary.txt chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt
      echo -e "No AMR genes found" > chroquetas_out/~{samplename}.ChroQueTaS.AMR_stats.txt
      echo "" | tee AMR_SUMMARY_STRING
      echo "" | tee ANNOTATED_AMR_SUMMARY_STRING
      echo "PASS" | tee CHROQUETAS_STATUS
    else
      echo "ERROR: Missing output files" | tee CHROQUETAS_STATUS
    fi
  fi
  >>>
  output {
    File? amr_stats_file = "chroquetas_out/~{samplename}.ChroQueTaS.AMR_stats.txt"
    File? amr_summary_file = "chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt"
    String chroquetas_mutations = read_string("AMR_SUMMARY_STRING")
    String chroquetas_fungicide_resistance = read_string("ANNOTATED_AMR_SUMMARY_STRING")
    String chroquetas_version = read_string("CHROQUETAS_VERSION")
    String chroquetas_status = read_string("CHROQUETAS_STATUS")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
  }
}