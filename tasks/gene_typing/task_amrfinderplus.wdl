version 1.0

task amrfinderplus_nuc {
  input {
    File assembly
    String samplename
    # Parameters 
    # --indent_min Minimum DNA %identity [0-1]; default is 0.9 (90%) or curated threshold if it exists
    # --mincov Minimum DNA %coverage [0-1]; default is 0.5 (50%)
    String? organism # make optional?
    Float? minid
    Float? mincov
    Boolean detailed_drug_class = false
    Int cpu = 4
    String docker = "staphb/ncbi-amrfinderplus:3.10.42"
    Int disk_size = 100
    Boolean hide_point_mutations = false
  }
  command <<<
    # logging info
    date | tee DATE
    amrfinder --version | tee AMRFINDER_VERSION
    
    # capture the database version; strip out unnecessary output, remove "Database version: " that prints in front of the actual database version
    amrfinder --database_version 2>/dev/null | grep "Database version" | sed 's|Database version: ||' | tee AMRFINDER_DB_VERSION

    ### set $amrfinder_organism BASH variable based on gambit_predicted_taxon or user-defined input string
    ### final variable has strict syntax/spelling based on list from amrfinder --list_organisms
    # there may be other Acinetobacter species to add later, like those in the A. baumannii-calcoaceticus species complex
    if [[ "~{organism}" == *"Acinetobacter"*"baumannii"* ]]; then
      amrfinder_organism="Acinetobacter_baumannii"
    elif [[ "~{organism}" == *"Campylobacter"*"coli"* ]] || [[ "~{organism}" == *"Campylobacter"*"jejuni"* ]]; then
      amrfinder_organism="Campylobacter"
    elif [[ "~{organism}" == *"Clostridioides"*"difficile"* ]]; then
      amrfinder_organism="Clostridioides_difficile"
    elif [[ "~{organism}" == *"Enterococcus"*"faecalis"* ]]; then 
      amrfinder_organism="Enterococcus_faecalis"
    elif [[ "~{organism}" == *"Enterococcus"*"faecium"* ]] || [[ "~{organism}" == *"Enterococcus"*"hirae"* ]]; then 
      amrfinder_organism="Enterococcus_faecium"
    # should capture all Shigella and Escherichia species
    elif [[ "~{organism}" == *"Escherichia"* ]] || [[ "~{organism}" == *"Shigella"* ]]; then 
      amrfinder_organism="Escherichia"
    # add other Klebsiella species later? Cannot use K. oxytoca as per amrfinderplus wiki
    elif [[ "~{organism}" == *"Klebsiella"*"aerogenes"* ]] || [[ "~{organism}" == *"Klebsiella"*"pnemoniae"* ]]; then 
      amrfinder_organism="Klebsiella"
    # because some people spell the species 'gonorrhea' differently
    elif [[ "~{organism}" == *"Neisseria"*"gonorrhea"* ]] || [[ "~{organism}" == *"Neisseria"*"gonorrhoeae"* ]] || [[ "~{organism}" == *"Neisseria"*"meningitidis"* ]]; then 
      amrfinder_organism="Neisseria"
    elif [[ "~{organism}" == *"Pseudomonas"*"aeruginosa"* ]]; then 
      amrfinder_organism="Pseudomonas_aeruginosa"
    # pretty broad, could work on Salmonella bongori and other species
    elif [[ "~{organism}" == *"Salmonella"* ]]; then 
      amrfinder_organism="Salmonella"
    elif [[ "~{organism}" == *"Staphylococcus"*"aureus"* ]]; then 
      amrfinder_organism="Staphylococcus_aureus"
    elif [[ "~{organism}" == *"Staphylococcus"*"pseudintermedius"* ]]; then 
      amrfinder_organism="Staphylococcus_pseudintermedius"
    elif [[ "~{organism}" == *"Streptococcus"*"agalactiae"* ]]; then 
      amrfinder_organism="Streptococcus_agalactiae"
    elif [[ "~{organism}" == *"Streptococcus"*"pneumoniae"* ]] || [[ "~{organism}" == *"Streptococcus"*"mitis"* ]]; then 
      amrfinder_organism="Streptococcus_pneumoniae"
    elif [[ "~{organism}" == *"Streptococcus"*"pyogenes"* ]]; then 
      amrfinder_organism="Streptococcus_pyogenes"
    elif [[ "~{organism}" == *"Vibrio"*"cholerae"* ]]; then 
      amrfinder_organism="Vibrio_cholerae"
    else 
      echo "Either Gambit predicted taxon is not supported by NCBI-AMRFinderPlus or the user did not supply an organism as input."
      echo "Skipping the use of amrfinder --organism optional parameter."
    fi

    # checking bash variable
    echo "amrfinder_organism is set to:" ${amrfinder_organism}
    
    # if amrfinder_organism variable is set, use --organism flag, otherwise do not use --organism flag
    if [[ -v amrfinder_organism ]] ; then
      # always use --plus flag, others may be left out if param is optional and not supplied 
      amrfinder --plus \
        --organism ${amrfinder_organism} \
        ~{'--name ' + samplename} \
        ~{'--nucleotide ' + assembly} \
        ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
        ~{'--threads ' + cpu} \
        ~{'--coverage_min ' + mincov} \
        ~{'--ident_min ' + minid}
    else 
      # always use --plus flag, others may be left out if param is optional and not supplied 
      amrfinder --plus \
        ~{'--name ' + samplename} \
        ~{'--nucleotide ' + assembly} \
        ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
        ~{'--threads ' + cpu} \
        ~{'--coverage_min ' + mincov} \
        ~{'--ident_min ' + minid}
    fi

    # remove mutations where Element subtype is "POINT"
    if [[ "~{hide_point_mutations}" == "true" ]]; then
      awk -F "\t" '$11 != "POINT"' ~{samplename}_amrfinder_all.tsv >> temp.tsv
      mv temp.tsv ~{samplename}_amrfinder_all.tsv
    fi

    # Element Type possibilities: AMR, STRESS, and VIRULENCE 
    # create headers for 3 output files; tee to 3 files and redirect STDOUT to dev null so it doesn't print to log file
    head -n 1 ~{samplename}_amrfinder_all.tsv | tee ~{samplename}_amrfinder_stress.tsv ~{samplename}_amrfinder_virulence.tsv ~{samplename}_amrfinder_amr.tsv >/dev/null
    # looks for all rows with STRESS, AMR, or VIRULENCE and append to TSVs
    grep 'STRESS' ~{samplename}_amrfinder_all.tsv >> ~{samplename}_amrfinder_stress.tsv
    grep 'VIRULENCE' ~{samplename}_amrfinder_all.tsv >> ~{samplename}_amrfinder_virulence.tsv
    # || true is so that the final grep exits with code 0, preventing failures
    grep 'AMR' ~{samplename}_amrfinder_all.tsv >> ~{samplename}_amrfinder_amr.tsv || true

    # create string outputs for all genes identified in AMR, STRESS, VIRULENCE
    amr_genes=$(awk -F '\t' '{ print $7 }' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//')
    stress_genes=$(awk -F '\t' '{ print $7 }' ~{samplename}_amrfinder_stress.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//')
    virulence_genes=$(awk -F '\t' '{ print $7 }' ~{samplename}_amrfinder_virulence.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//')
    
    if [[ "~{detailed_drug_class}" == "true" ]]; then
      # create string outputs for AMR drug classes
      amr_classes=$(awk -F '\t' 'BEGIN{OFS=":"} {print $7,$12}' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//')
      # create string outputs for AMR drug subclasses
      amr_subclasses=$(awk -F '\t' 'BEGIN{OFS=":"} {print $7,$13}' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//')
    else
      amr_classes=$(awk -F '\t' '{ print $12 }' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | uniq | tr '\n' ', ' | sed 's/.$//')
      amr_subclasses=$(awk -F '\t' '{ print $13 }' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | uniq | tr '\n' ', ' | sed 's/.$//')
    fi

    # if variable for list of genes is EMPTY, write string saying it is empty to float to Terra table
    if [ -z "${amr_genes}" ]; then
       amr_genes="No AMR genes detected by NCBI-AMRFinderPlus"
    fi 
    if [ -z "${stress_genes}" ]; then
       stress_genes="No STRESS genes detected by NCBI-AMRFinderPlus"
    fi 
    if [ -z "${virulence_genes}" ]; then
       virulence_genes="No VIRULENCE genes detected by NCBI-AMRFinderPlus"
    fi 
    if [ -z "${amr_classes}" ]; then
       amr_classes="No AMR genes detected by NCBI-AMRFinderPlus"
    fi 
    if [ -z "${amr_subclasses}" ]; then
       amr_subclasses="No AMR genes detected by NCBI-AMRFinderPlus"
    fi 

    # create final output strings
    echo "${amr_genes}" > AMR_GENES
    echo "${stress_genes}" > STRESS_GENES
    echo "${virulence_genes}" > VIRULENCE_GENES
    echo "${amr_classes}" > AMR_CLASSES
    echo "${amr_subclasses}" > AMR_SUBCLASSES
  >>>
  output {
    File amrfinderplus_all_report = "~{samplename}_amrfinder_all.tsv"
    File amrfinderplus_amr_report = "~{samplename}_amrfinder_amr.tsv"
    File amrfinderplus_stress_report = "~{samplename}_amrfinder_stress.tsv"
    File amrfinderplus_virulence_report = "~{samplename}_amrfinder_virulence.tsv"
    String amrfinderplus_amr_genes = read_string("AMR_GENES")
    String amrfinderplus_stress_genes = read_string("STRESS_GENES")
    String amrfinderplus_virulence_genes = read_string("VIRULENCE_GENES")
    String amrfinderplus_amr_classes = read_string("AMR_CLASSES")
    String amrfinderplus_amr_subclasses = read_string("AMR_SUBCLASSES")
    String amrfinderplus_version = read_string("AMRFINDER_VERSION")
    String amrfinderplus_db_version = read_string("AMRFINDER_DB_VERSION")
  }
  runtime {
    memory: "8 GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
