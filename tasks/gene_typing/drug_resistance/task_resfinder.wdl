version 1.0

task resfinder {
  input {
    File assembly # Input fasta file
    String samplename
    String? organism # Species in the sample, species should be entered with their full scientific names (e.g. "escherichia coli"), using quotation marks
    Boolean acquired = true # Run resfinder for acquired resistance genes
    Float min_cov = 0.5 # Minimum (breadth-of) coverage of ResFinder
    Float min_id = 0.9 # Threshold for identity of ResFinder
    Boolean call_pointfinder = false # Run pointfinder for chromosomal mutations
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/resfinder:4.1.11"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    date | tee DATE
    run_resfinder.py --version | tee RESFINDER_VERSION
    echo "unmodified from resfinder docker container" > RESFINDER_DB_VERSION

    # set $resfinder_organism BASH variable based on gambit_predicted_taxon or user-defined input string
    if [[ "~{organism}" == *"Campylobacter"*"jejuni"* ]]; then
      resfinder_organism="campylobacter jejuni"
    elif [[ "~{organism}" == *"Campylobacter"*"coli"* ]]; then
      resfinder_organism="campylobacter coli"
    elif [[ "~{organism}" == *"Campylobacter"* ]]; then
      resfinder_organism="campylobacter"
    elif [[ "~{organism}" == *"Enterococcus"*"faecalis"* ]]; then 
      resfinder_organism="enterococcus faecalis"
    elif [[ "~{organism}" == *"Enterococcus"*"faecium"* ]]; then 
      resfinder_organism="enterococcus faecium"
    # to allow for both E. coli and any Shigella species to be processed via PointFinder as E. coli
    elif [[ "~{organism}" == *"Escherichia"*"coli"* ]] || [[ "~{organism}" == *"Shigella"* ]]; then 
      resfinder_organism="escherichia coli"
    elif [[ "~{organism}" == *"Klebsiella"* ]]; then 
      resfinder_organism="klebsiella"
    elif [[ "~{organism}" == *"Neisseria"*"gonorrhoeae"* ]]; then 
      resfinder_organism="neisseria gonorrhoeae"
    elif [[ "~{organism}" == *"Salmonella"* ]]; then 
      resfinder_organism="salmonella"
    elif [[ "~{organism}" == *"Staphylococcus"*"aureus"* ]]; then 
      resfinder_organism="staphylococcus aureus"
    elif [[ "~{organism}" == *"Mycobacterium"*"tuberculosis"* ]]; then 
      resfinder_organism="mycobacterium tuberculosis"
    elif [[ "~{organism}" == *"Helicobacter"*"pylori"* ]]; then 
      resfinder_organism="helicobacter pylori"
    else 
      echo "Either Gambit predicted taxon is not supported by resfinder or the user did not supply an organism as input."
      echo "Skipping the use of resfinder --species optional parameter."
      echo "WARNING: This will disable PointFinder due to the requirement of --species flag."
    fi

    # if resfinder_organism variable is set, use --species flag, otherwise do not use --species flag
    if [[ -v resfinder_organism ]] ; then
      run_resfinder.py \
        --inputfasta ~{assembly} \
        --outputPath . \
        --species "${resfinder_organism}" \
        ~{true="--acquired" false="" acquired} \
        ~{'--min_cov ' + min_cov} \
        ~{'--threshold ' + min_id} \
        ~{true="--point" false="" call_pointfinder} 
    else 
      # pointfinder requires the use of the --species flag, so if resfinder_organism is not set, do not run pointfinder
      run_resfinder.py \
        --inputfasta ~{assembly} \
        --outputPath . \
        --species "other" \
        ~{true="--acquired" false="" acquired} \
        ~{'--min_cov ' + min_cov} \
        ~{'--threshold ' + min_id}
    fi

    # replace space in resfinder_organism with underscore
    resfinder_organism="${resfinder_organism// /_}"

    # rename files
    mv -v pheno_table.txt ~{samplename}_pheno_table.tsv
    if [ -f "pheno_table_${resfinder_organism}.txt" ]; then
      # rename file to have proper extension & samplename included
      mv -v "pheno_table_${resfinder_organism}.txt" ~{samplename}_pheno_table_species.tsv
    fi
    mv -v ResFinder_Hit_in_genome_seq.fsa ~{samplename}_ResFinder_Hit_in_genome_seq.fsa
    mv -v ResFinder_Resistance_gene_seq.fsa ~{samplename}_ResFinder_Resistance_gene_seq.fsa
    mv -v ResFinder_results_tab.txt ~{samplename}_ResFinder_results_tab.tsv

    # if pointfinder was run, rename files
    if [ -f PointFinder_prediction.txt ]; then
      mv -v PointFinder_prediction.txt ~{samplename}_PointFinder_prediction.tsv
      mv -v PointFinder_results.txt ~{samplename}_PointFinder_results.tsv
    fi

    # parse ~{samplename}_pheno_table.tsv for predicted phenotypes and genes associated
    # strip off 18 lines from top of file (18th line is the header with the columns: antimicrobial, class, WGS-predicted phenotype, Match, Genetic Background)
    tail +18 ~{samplename}_pheno_table.tsv > ~{samplename}_pheno_table.headerless.tsv
    
    # convert all letters in first column (antibiotic) to uppercase. For readability of output string
    awk -F '\t' 'BEGIN{OFS="\t"} { $1=toupper($1) } 1' ~{samplename}_pheno_table.headerless.tsv > ~{samplename}_pheno_table.headerless.uppercase.tsv

    # if column 3 shows 'Resistant', then print list of drugs followed by the genes/point mutations responsible
    awk -F '\t' 'BEGIN{OFS=":"; ORS="; "} { if($3 == "Resistant") {print $1,$5}}' ~{samplename}_pheno_table.headerless.uppercase.tsv \
    | sed 's/..$//' > RESFINDER_PREDICTED_PHENO_RESISTANCE.txt

    # check for XDR Shigella status, based on CDC definition here: https://emergency.cdc.gov/han/2023/han00486.asp
    # requirements:
    # organism input (i.e. gambit_predicted_taxon) must contain "Shigella"
    # predicted resistance to antimicrobials must include ALL: "ceftriaxone", "azithromycin", "ciprofloxacin", "trimethoprim", "sulfamethoxazole", and "ampicillin"
    if [[ "~{organism}" != *"Shigella"* ]]; then
      echo 'Either the input gambit_predicted_taxon does not contain "Shigella" or the user did not supply the organism as an input string to the workflow.'
      echo "Skipping XDR Shigella check."
      echo "Not Shigella based on gambit_predicted_taxon or user input" | tee RESFINDER_PREDICTED_XDR_SHIGELLA.txt
    # if organism input string DOES contain the word "Shigella", check for resistance predictions to 6 drugs in XDR definition
    elif [[ "~{organism}" == *"Shigella"* ]]; then
      # nested if: if grep finds the all drugs, set output to XDR shigella, but if not, set it to "Not XDR Shigella"
      if grep -qi "ceftriaxone" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt && \
      grep -qi "azithromycin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt && \
      grep -qi "ciprofloxacin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt && \
      grep -qi "trimethoprim" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt && \
      grep -qi "sulfamethoxazole" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt && \
      grep -qi "ampicillin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
        echo "XDR Shigella based on predicted resistance to ceftriaxone, azithromycin, ciprofloxacin, trimethoprim, sulfamethoxazole, and ampicillin. Please verify by reviewing ~{samplename}_pheno_table.tsv and ~{samplename}_ResFinder_results_tab.tsv"
        echo "XDR Shigella" > RESFINDER_PREDICTED_XDR_SHIGELLA.txt
      # ~{organism} does contain the word "Shigella", but one of the greps failed, meaning not all drug resistances' were predicted
      else
        echo "Not XDR Shigella" | tee RESFINDER_PREDICTED_XDR_SHIGELLA.txt
      fi
    fi

    # set output strings for Resistance or no Resistance predicted for each of 6 drugs
    # if grep finds the drug in the RESFINDER_PREDICTED_PHENO_RESISTANCE.txt file, then set the output string to "Resistance"
    # if grep does not find the drug in the RESFINDER_PREDICTED_PHENO_RESISTANCE.txt file, then set the output string to "No resistance predicted"
    # ampicillin
    if grep -qi "ampicillin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_AMP.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_AMP.txt
    fi
    # azithromycin
    if grep -qi "azithromycin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_AZM.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_AZM.txt
    fi
    # ceftriaxone
    if grep -qi "ceftriaxone" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_AXO.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_AXO.txt
    fi
    # ciprofloxacin
    if grep -qi "ciprofloxacin" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_CIP.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_CIP.txt
    fi
    # sulfamethoxazole
    if grep -qi "sulfamethoxazole" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_SMX.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_SMX.txt
    fi
    # trimethoprim
    if grep -qi "trimethoprim" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
      echo "Resistance" > RESFINDER_PREDICTED_RESISTANCE_TMP.txt
    else
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_TMP.txt
    fi

  >>>
  output {
    File resfinder_pheno_table = "~{samplename}_pheno_table.tsv"
    File? resfinder_pheno_table_species = "~{samplename}_pheno_table_species.tsv"
    File resfinder_hit_in_genome_seq = "~{samplename}_ResFinder_Hit_in_genome_seq.fsa"
    File resfinder_resistance_gene_seq = "~{samplename}_ResFinder_Resistance_gene_seq.fsa"
    File resfinder_results_tab = "~{samplename}_ResFinder_results_tab.tsv"
    File? pointfinder_pheno_table = "~{samplename}_PointFinder_prediction.tsv"
    File? pointfinder_results = "~{samplename}_PointFinder_results.tsv"
    String resfinder_predicted_pheno_resistance = read_string("RESFINDER_PREDICTED_PHENO_RESISTANCE.txt")
    String resfinder_predicted_xdr_shigella = read_string("RESFINDER_PREDICTED_XDR_SHIGELLA.txt")
    String resfinder_predicted_resistance_Amp = read_string("RESFINDER_PREDICTED_RESISTANCE_AMP.txt")
    String resfinder_predicted_resistance_Azm = read_string("RESFINDER_PREDICTED_RESISTANCE_AZM.txt")
    String resfinder_predicted_resistance_Axo = read_string("RESFINDER_PREDICTED_RESISTANCE_AXO.txt")
    String resfinder_predicted_resistance_Cip = read_string("RESFINDER_PREDICTED_RESISTANCE_CIP.txt")
    String resfinder_predicted_resistance_Smx = read_string("RESFINDER_PREDICTED_RESISTANCE_SMX.txt")
    String resfinder_predicted_resistance_Tmp = read_string("RESFINDER_PREDICTED_RESISTANCE_TMP.txt")
    String resfinder_docker = "~{docker}"
    String resfinder_version = read_string("RESFINDER_VERSION")
    String resfinder_db_version = read_string("RESFINDER_DB_VERSION")
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3  
  }
}
