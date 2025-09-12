version 1.0

task chroquetas {
  input {
    File genome_fasta
    String species
    String samplename

    Int min_cov = 75 # set to ChroQueTas default
    Int min_id = 90  # set to ChroQueTas default
    File? proteome_fasta
    String? translation_code

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/chroquetas:1.0.0"
    Int cpu = 2
    Int disk_size = 16
    Int memory = 8
  }
  command <<<
  # get version (non-zero exit status)
  ChroQueTas.sh --version | sed -E 's/\(Chromosome Query Targets\) version //' | tee CHROQUETAS_VERSION

  # fail hard
  set -euo pipefail

  # ChroQueTas expects species with underscore delimiters and a capital first letter
  species_prep=$(echo ~{species} | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
  corrected_species=$(echo ${species_prep^})

  # call chroquetas
  ChroQueTas.sh \
    -g ~{genome_fasta} \
    -s ${corrected_species} \
    --min_cov ~{min_cov} \
    --min_id ~{min_id} \
    -t ~{cpu} \
    -o chroquetas_out \
    ~{if defined(proteome_fasta) then "-p ~{proteome_fasta}" else ""} \
    ~{if defined(translation_code) then "-c ~{translation_code}" else ""}

  # rename output files to include sample name
  based_name=$(echo $(basename ~{genome_fasta}) | sed -E 's/(.*)\.[^\.]+$/\1/')
  mv chroquetas_out/${based_name}.ChroQueTaS.AMR_stats.txt chroquetas_out/~{samplename}.ChroQueTaS.AMR_stats.txt
  mv chroquetas_out/${based_name}.ChroQueTaS.AMR_summary.txt chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt

  # extract AMR summary string
  tail -n+2 chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt \
    | awk '{ print $1, $4, $3, $5 }' \
    | sed -E 's/([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)/\1_\2\3\4/' \
    | tr '\n' ',' \
    | sed -E 's/,$//' \
    | tee AMR_SUMMARY_STRING
  >>>
  output {
    File amr_stats_file = "chroquetas_out/~{samplename}.ChroQueTaS.AMR_stats.txt"
    File amr_summary_file = "chroquetas_out/~{samplename}.ChroQueTaS.AMR_summary.txt"
    String amr_summary_string = read_string("AMR_SUMMARY_STRING")
    String chroquetas_version = read_string("CHROQUETAS_VERSION")
  }

  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
  }
}