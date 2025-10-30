version 1.0

task amrfinderplus_nuc {
  input {
    File assembly
    File? annotation_assembly
    File? protein_fasta
    File? gff
    Boolean? use_gff = false
    String samplename
    # Parameters 
    # --indent_min Minimum DNA %identity [0-1]; default is 0.9 (90%) or curated threshold if it exists
    # --mincov Minimum DNA %coverage [0-1]; default is 0.5 (50%)
    String? organism 
    String? annotation_format
    Float? min_percent_identity
    Float? min_percent_coverage
    Boolean detailed_drug_class = false
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-amrfinderplus:4.0.23-2025-07-16.1"
    Int disk_size = 50
    Int memory = 8
    Boolean hide_point_mutations = false
    Boolean separate_betalactam_genes = false
  }
  command <<<
    # logging info
    date | tee DATE
    amrfinder --version | tee AMRFINDER_VERSION
    
    # capture the database version; strip out unnecessary output, remove "Database version: " that prints in front of the actual database version
    amrfinder --database_version 2>/dev/null | grep "Database version" | sed 's|Database version: ||' | tee AMRFINDER_DB_VERSION

    ## create associative array
    declare -A organisms=( ["Acinetobacter baumannii"]="Acinetobacter_baumannii" ["Burkholderia cepacia"]="Burkholderia_cepacia" \
      ["Burkholderia pseudomallei"]="Burkholderia_pseudomallei" ["Campylobacter coli"]="Campylobacter" ["Campylobacter jejuni"]="Campylobacter" \
      ["Citrobacter freundii"]="Citrobacter_freundii" ["Clostridioides difficile"]="Clostridioides_difficile" ["Enterobacter asburiae"]="Enterobacter_asburiae" ["Enterobacter cloacae"]="Enterobacter_cloacae" \
      ["Enterococcus faecalis"]="Enterococcus_faecalis" ["Enterococcus hirae"]="Enterococcus_faecium" ["Enterococcusfaecium"]="Enterococcus_faecium" \
      [*"Escherichia"*]="Escherichia" [*"Shigella"*]="Escherichia" ["Klebsiella aerogenes"]="Klebsiella_pneumoniae" ["Klebsiella pneumoniae"]="Klebsiella_pneumoniae" \
      ["Klebsiella variicola"]="Klebsiella_pneumoniae" ["Klebsiella oxytoca"]="Klebsiella_oxytoca" ["Neisseria gonorrhea"]="Neisseria_gonorrhoeae" \
      ["Neisseria gonorrhoeae"]="Neisseria_gonorrhoeae" ["Neisseria meningitidis"]="Neisseria_meningitidis" ["Pseudomonas aeruginosa"]="Pseudomonas_aeruginosa" \
      [*"Salmonella"*]="Salmonella" ["Serratia marcescens"]="Serratia_marcescens" ["Staphylococcus aureus"]="Staphylococcus_aureus" \
      ["Staphylococcus pseudintermedius"]="Staphylococcus_pseudintermedius" ["Streptococcus agalactiae"]="Streptococcus_agalactiae" \
      ["Streptococcus pneumoniae"]="Streptococcus_pneumoniae" ["Streptococcus mitis"]="Streptococcus_pneumoniae" \
      ["Streptococcus pyogenes"]="Streptococcus_pyogenes" ["Vibrio cholerae"]="Vibrio_cholerae" ["Vibrio parahaemolyticus"]="Vibrio_parahaemolyticus" ["Vibrio vulnificus"]="Vibrio_vulnificus"
      )

    for key in "${!organisms[@]}"; do
      if [[ "~{organism}" == $key ]]; then
        amrfinder_organism=${organisms[$key]}
      elif [[ "~{organism}" == *$key* ]]; then
        amrfinder_organism=${organisms[$key]}
      fi
    done

    # checking bash variable
    echo "amrfinder_organism is set to: ${amrfinder_organism}"
    
    # if amrfinder_organism variable is set, use --organism flag, otherwise do not use --organism flag
    if [[ -v amrfinder_organism ]] ; then
      # always use --plus flag, others may be left out if param is optional and not supplied
      if [[ "~{use_gff}" == "true" ]]; then
        amrfinder --plus \
          --organism "${amrfinder_organism}" \
          ~{'--name ' + samplename} \
          ~{'--nucleotide ' + annotation_assembly} \
          ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
          ~{'--threads ' + cpu} \
          ~{'--protein ' + protein_fasta} \
          ~{'--gff ' + gff} \
          ~{'--annotation_format '+ annotation_format} \
          ~{'--coverage_min ' + min_percent_coverage} \
          ~{'--ident_min ' + min_percent_identity}
      else
        amrfinder --plus \
          --organism "${amrfinder_organism}" \
          ~{'--name ' + samplename} \
          ~{'--nucleotide ' + assembly} \
          ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
          ~{'--threads ' + cpu} \
          ~{'--coverage_min ' + min_percent_coverage} \
          ~{'--ident_min ' + min_percent_identity}
      fi
    else 
      echo "Either the organism (~{organism}) is not recognized by NCBI-AMRFinderPlus or the user did not supply an organism as input."
      echo "Skipping the use of amrfinder --organism optional parameter."
      # always use --plus flag, others may be left out if param is optional and not supplied
      if [[ "~{use_gff}" == "true" ]]; then
        amrfinder --plus \
          ~{'--name ' + samplename} \
          ~{'--nucleotide ' + annotation_assembly} \
          ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
          ~{'--threads ' + cpu} \
          ~{'--protein ' + protein_fasta} \
          ~{'--gff ' + gff} \
          ~{'--annotation_format '+ annotation_format} \
          ~{'--coverage_min ' + min_percent_coverage} \
          ~{'--ident_min ' + min_percent_identity}
      else
        amrfinder --plus \
          ~{'--name ' + samplename} \
          ~{'--nucleotide ' + assembly} \
          ~{'-o ' + samplename + '_amrfinder_all.tsv'} \
          ~{'--threads ' + cpu} \
          ~{'--coverage_min ' + min_percent_coverage} \
          ~{'--ident_min ' + min_percent_identity}
      fi      
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
    amr_core_genes=$(awk -F '\t' '{ if($9 == "core") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//')
    amr_plus_genes=$(awk -F '\t' '{ if($9 != "core") { print $7}}' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | tr '\n' ', ' | sed 's/.$//')
    stress_genes=$(awk -F '\t' '{ print $7 }' ~{samplename}_amrfinder_stress.tsv | tail -n+2 | sort | tr '\n' ', ' | sed 's/.$//')
    virulence_genes=$(awk -F '\t' '{ print $7 }' ~{samplename}_amrfinder_virulence.tsv | tail -n+2 | sort | tr '\n' ', ' | sed 's/.$//')

    if [[ "~{detailed_drug_class}" == "true" ]]; then
      # create string outputs for AMR drug classes
      amr_classes=$(awk -F '\t' 'BEGIN{OFS=":"} {print $7,$12}' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | tr '\n' ', ' | sed 's/.$//')
      # create string outputs for AMR drug subclasses
      amr_subclasses=$(awk -F '\t' 'BEGIN{OFS=":"} {print $7,$13}' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | tr '\n' ', ' | sed 's/.$//')
    else
      amr_classes=$(awk -F '\t' '{ print $12 }' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | uniq | tr '\n' ', ' | sed 's/.$//')
      amr_subclasses=$(awk -F '\t' '{ print $13 }' ~{samplename}_amrfinder_amr.tsv | tail -n+2 | sort | uniq | tr '\n' ', ' | sed 's/.$//')
    fi

    if [[ "~{separate_betalactam_genes}" == "true" ]]; then
      betalactam_genes=$(awk -F '\t' '{ if($12 == "BETA-LACTAM") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//') 
      betalactam_betalactam_genes=$(awk -F '\t' '{ if($13 == "BETA-LACTAM") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//') 
      betalactam_carbapenem_genes=$(awk -F '\t' '{ if($13 == "CARBAPENEM") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//')
      betalactam_cephalosporin_genes=$(awk -F '\t' '{ if($13 == "CEPHALOSPORIN") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//')
      betalactam_cephalothin_genes=$(awk -F '\t' '{ if($13 == "CEPHALOTHIN") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//')
      betalactam_methicillin_genes=$(awk -F '\t' '{ if($13 == "METHICILLIN") { print $7}}' ~{samplename}_amrfinder_amr.tsv | sort | tr '\n' ', ' | sed 's/.$//')

      # if variable for list of genes is EMPTY, write string saying it is empty to float to Terra table
      if [ -z "${betalactam_genes}" ]; then
        betalactam_genes="No BETA-LACTAM genes detected by NCBI-AMRFinderPlus"
      fi
      if [ -z "${betalactam_betalactam_genes}" ]; then
        betalactam_betalactam_genes="No BETA-LACTAM genes detected by NCBI-AMRFinderPlus"
      fi
      if [ -z "${betalactam_carbapenem_genes}" ]; then
        betalactam_carbapenem_genes="No BETA-LACTAM CARBAPENEM genes detected by NCBI-AMRFinderPlus"
      fi
      if [ -z "${betalactam_cephalosporin_genes}" ]; then
        betalactam_cephalosporin_genes="No BETA-LACTAM CEPHALOSPORIN genes detected by NCBI-AMRFinderPlus"
      fi
      if [ -z "${betalactam_cephalothin_genes}" ]; then
        betalactam_cephalothin_genes="No BETA-LACTAM CEPHALOTHIN genes detected by NCBI-AMRFinderPlus"
      fi
      if [ -z "${betalactam_methicillin_genes}" ]; then
        betalactam_methicillin_genes="No BETA-LACTAM METHICILLIN genes detected by NCBI-AMRFinderPlus"
      fi
    fi

    # if variable for list of genes is EMPTY, write string saying it is empty to float to Terra table
    if [ -z "${amr_core_genes}" ]; then
       amr_core_genes="No core AMR genes detected by NCBI-AMRFinderPlus"
    fi 
    if [ -z "${amr_plus_genes}" ]; then
       amr_plus_genes="No plus AMR genes detected by NCBI-AMRFinderPlus"
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
    echo "${amr_core_genes}" > AMR_CORE_GENES
    echo "${amr_plus_genes}" > AMR_PLUS_GENES
    echo "${stress_genes}" > STRESS_GENES
    echo "${virulence_genes}" > VIRULENCE_GENES
    echo "${amr_classes}" > AMR_CLASSES
    echo "${amr_subclasses}" > AMR_SUBCLASSES
  
    # if separate_betalactam_genes is true, create final output strings (if false, these values will be blank)
    echo "${betalactam_genes}" > BETA_LACTAM_GENES
    echo "${betalactam_betalactam_genes}" > BETA_LACTAM_BETA_LACTAM_GENES
    echo "${betalactam_carbapenem_genes}" > BETA_LACTAM_CARBAPENEM_GENES
    echo "${betalactam_cephalosporin_genes}" > BETA_LACTAM_CEPHALOSPORIN_GENES
    echo "${betalactam_cephalothin_genes}" > BETA_LACTAM_CEPHALOTHIN_GENES
    echo "${betalactam_methicillin_genes}" > BETA_LACTAM_METHICILLIN_GENES
  >>>
  output {
    File amrfinderplus_all_report = "~{samplename}_amrfinder_all.tsv"
    File amrfinderplus_amr_report = "~{samplename}_amrfinder_amr.tsv"
    File amrfinderplus_stress_report = "~{samplename}_amrfinder_stress.tsv"
    File amrfinderplus_virulence_report = "~{samplename}_amrfinder_virulence.tsv"
    String amrfinderplus_amr_core_genes = read_string("AMR_CORE_GENES")
    String amrfinderplus_amr_plus_genes = read_string("AMR_PLUS_GENES")
    String amrfinderplus_stress_genes = read_string("STRESS_GENES")
    String amrfinderplus_virulence_genes = read_string("VIRULENCE_GENES")
    String amrfinderplus_amr_classes = read_string("AMR_CLASSES")
    String amrfinderplus_amr_subclasses = read_string("AMR_SUBCLASSES")
    String amrfinderplus_version = read_string("AMRFINDER_VERSION")
    String amrfinderplus_db_version = read_string("AMRFINDER_DB_VERSION")

    String amrfinderplus_amr_betalactam_genes = read_string("BETA_LACTAM_GENES")
    String amrfinderplus_amr_betalactam_betalactam_genes = read_string("BETA_LACTAM_BETA_LACTAM_GENES")
    String amrfinderplus_amr_betalactam_carbapenem_genes = read_string("BETA_LACTAM_CARBAPENEM_GENES")
    String amrfinderplus_amr_betalactam_cephalosporin_genes = read_string("BETA_LACTAM_CEPHALOSPORIN_GENES")
    String amrfinderplus_amr_betalactam_cephalothin_genes = read_string("BETA_LACTAM_CEPHALOTHIN_GENES")
    String amrfinderplus_amr_betalactam_methicillin_genes = read_string("BETA_LACTAM_METHICILLIN_GENES")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}
