version 1.0

task midas {
  input {
    File read1
    File? read2
    File midas_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz"
    String samplename
    String docker = "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    Int? memory = 32
    Int? cpu = 4
  }
  command <<<
    date | tee DATE

    # Decompress the Midas database
    mkdir db
    tar -C ./db/ -xzvf ~{midas_db}  

    # Run Midas
    run_midas.py species ~{samplename} -1 ~{read1} ~{'-2 ' + read2} -d db/midas_db_v1.2/ -t ~{cpu} 

    # rename output files
    mv ~{samplename}/species/species_profile.txt ~{samplename}/species/~{samplename}_species_profile.tsv
    mv ~{samplename}/species/log.txt ~{samplename}/species/~{samplename}_log.txt

    # determine if secondary species
    # remove header
    cat ~{samplename}/species/~{samplename}_species_profile.tsv | tail -n +2 > output.tsv
    # get primary genus: sort by coverage (descending), get top non-header row, cut for species_ID column, parse column to get only genus name
    primary_genus=$(cat output.tsv | sort -k 3 -r -n | awk 'NR==1' | cut -f1 | cut -f1 -d"_")

    # filter to remove lines with primary genus
    grep -v -F "$primary_genus" output.tsv > output2.tsv

    # get secondary species: sort by coverage again to be safe, get top non-header row, cut for species_ID column, parse column to get only genus name
    secondary_genus=$(cat output2.tsv | sort -k 3 -r -n | awk 'NR==1' | cut -f1 | cut -f1 -d"_")
    # get coverage of secondary genus
    secondary_genus_coverage=$(cat output2.tsv | sort -k 3 -r -n | awk 'NR==1' | cut -f3 )
    # round coverage of secondary genus to three decimal places
    secondary_genus_coverage=$(printf %.3f $secondary_genus_coverage)

    # indicate if no secondary genus was detected
    cov_less_than_one=$(echo ${secondary_genus_coverage} 1.0 | awk '{if ($1 < $2) print "true"; else print "false"}')
    if [[ "${cov_less_than_one}" == "true" ]] ; then
       secondary_genus="No secondary genus detected (>1.0X coverage)"
    fi
    
    # create final output strings
    echo "${primary_genus}" > PRIMARY_GENUS
    echo "${secondary_genus}" > SECONDARY_GENUS
    echo "${secondary_genus_coverage}" > SECONDARY_GENUS_COVERAGE

  >>>
  output {
    String midas_docker = docker
    String midas_analysis_date = read_string("DATE")
    File midas_report = "~{samplename}/species/~{samplename}_species_profile.tsv"
    File midas_log = "~{samplename}/species/~{samplename}_log.txt"
    String midas_primary_genus = read_string("PRIMARY_GENUS")
    String midas_secondary_genus = read_string("SECONDARY_GENUS")
    String midas_secondary_genus_coverage = read_string("SECONDARY_GENUS_COVERAGE")
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}