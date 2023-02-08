version 1.0

task midas_theiaeuk {
  input {
    File read1
    File? read2
    File midas_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz"
    String samplename
    String docker = "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    Int memory = 32
    Int cpu = 4
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

task midas_theiaprok {
  input {
    File read1
    File? read2
    File midas_db = "gs://theiagen-public-files-rp/terra/theiaprok-files/midas/midas_db_v1.2.tar.gz"
    Int disk_size = 100
    String samplename
    String docker = "quay.io/fhcrc-microbiome/midas:v1.3.2--6"
    Int memory = 32
    Int cpu = 4
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

    # Run a python block to parse output file for terra data tables
    # pandas is available in default docker image for python2 but not python3
    python2 <<CODE
    import pandas as pd

    df = pd.read_csv('~{samplename}/species/~{samplename}_species_profile.tsv',sep="\t", header=0)
    # round relative abundance to 4 decimal places
    df = df.round(4)
    # sort by relative abundance
    sorted_df = df.sort_values(by=['relative_abundance'], ascending=False)
    # split species_id column 
    sorted_df[['genus','species','strain']] = sorted_df['species_id'].str.split(pat = '_',expand=True,n=2)
    # capture primary genus
    primary_genus = sorted_df.genus.iloc[0]
    # remove rows where genus is primary_genus 
    filtered_df = sorted_df[sorted_df['genus'].str.contains(str(primary_genus))==False ]
    # re-sort by relative abundance, just in case
    filtered_sorted_df = filtered_df.sort_values(by=['relative_abundance'], ascending=False)
    # capture secondary genus
    secondary_genus = filtered_sorted_df.genus.iloc[0]
    # capture abundance of secondary genus
    secondary_genus_abundance = filtered_sorted_df.relative_abundance.iloc[0]
    # if secondary genus abundance is less than one, replace genus with text indicating no secondary genus detected
    if secondary_genus_abundance < 0.01:
      secondary_genus="No secondary genus detected (>1% relative abundance)"

    # write text files
    with open("PRIMARY_GENUS", 'wt') as pg:
      pg.write(str(primary_genus))
    with open("SECONDARY_GENUS", 'wt') as sg:
      sg.write(str(secondary_genus))
    with open("SECONDARY_GENUS_ABUNDANCE", 'wt') as sga:
      sga.write(str(secondary_genus_abundance))

    CODE
  >>>
  output {
    String midas_docker = docker
    String midas_analysis_date = read_string("DATE")
    File midas_report = "~{samplename}/species/~{samplename}_species_profile.tsv"
    File midas_log = "~{samplename}/species/~{samplename}_log.txt"
    String midas_primary_genus = read_string("PRIMARY_GENUS")
    String midas_secondary_genus = read_string("SECONDARY_GENUS")
    Float midas_secondary_genus_abundance = read_string("SECONDARY_GENUS_ABUNDANCE")
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