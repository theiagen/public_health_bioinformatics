version 1.0

task resfinder {
  input {
    File assembly
    String samplename
    String? organism # Species in the sample, species should be entered with their full scientific names (e.g. "escherichia coli"), using quotation marks
    Boolean acquired = true # Run resfinder for acquired resistance genes
    Float min_percent_coverage = 0.5 # Minimum (breadth-of) coverage of ResFinder
    Float min_percent_identity = 0.9 # Threshold for identity of ResFinder
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

    # set resfinder_organism using an associative array and loop
    # more specific keys should come before general keys to avoid premature matches
    declare -A organism_map=(
      ["Campylobacter jejuni"]="campylobacter jejuni"
      ["Campylobacter coli"]="campylobacter coli"
      ["Campylobacter"]="campylobacter"
      ["Enterococcus faecalis"]="enterococcus faecalis"
      ["Enterococcus faecium"]="enterococcus faecium"
      ["Escherichia coli"]="escherichia coli"
      ["Shigella"]="escherichia coli"
      ["Klebsiella"]="klebsiella"
      ["Neisseria gonorrhoeae"]="neisseria gonorrhoeae"
      ["Salmonella"]="salmonella"
      ["Staphylococcus aureus"]="staphylococcus aureus"
      ["Mycobacterium tuberculosis"]="mycobacterium tuberculosis"
      ["Helicobacter pylori"]="helicobacter pylori"
    )

    resfinder_organism=""
    for key in "${!organism_map[@]}"; do
      if [[ "~{organism}" =~ $key ]]; then
        resfinder_organism="${organism_map[$key]}"
        break
      fi
    done

    # run resfinder with either resfinder_organism and pointfinder, or not 
    if [[ -z "$resfinder_organism" ]]; then
      echo "Either Gambit predicted taxon is not supported by resfinder or the user did not supply an organism as input."
      echo "Skipping the use of resfinder --species optional parameter."
      echo "WARNING: This will disable PointFinder due to the requirement of --species flag."
  
      run_resfinder.py \
        --inputfasta ~{assembly} \
        --outputPath . \
        --species "other" \
        ~{true="--acquired" false="" acquired} \
        ~{'--min_cov ' + min_percent_coverage} \
        ~{'--threshold ' + min_percent_identity}
    else 

      run_resfinder.py \
        --inputfasta ~{assembly} \
        --outputPath . \
        --species "${resfinder_organism}" \
        ~{true="--acquired" false="" acquired} \
        ~{'--min_cov ' + min_percent_coverage} \
        ~{'--threshold ' + min_percent_identity} \
        ~{true="--point" false="" call_pointfinder}

    fi

    # rename all output files using an associative array
    declare -A file_rename_map=(
      ["pheno_table.txt"]="~{samplename}_pheno_table.tsv"
      ["pheno_table_${resfinder_organism// /_}.txt"]="~{samplename}_pheno_table_species.tsv"
      ["ResFinder_Hit_in_genome_seq.fsa"]="~{samplename}_ResFinder_Hit_in_genome_seq.fsa"
      ["ResFinder_Resistance_gene_seq.fsa"]="~{samplename}_ResFinder_Resistance_gene_seq.fsa"
      ["ResFinder_results_tab.txt"]="~{samplename}_ResFinder_results_tab.tsv"
      ["PointFinder_prediction.txt"]="~{samplename}_PointFinder_prediction.tsv"
      ["PointFinder_results.txt"]="~{samplename}_PointFinder_results.tsv"
    )

    # check if file exists before renaming
    for file in "${!file_rename_map[@]}"; do
      if [ -f "$file" ]; then
        mv -v "$file" "${file_rename_map[$file]}"
      fi
    done

    # create an uppercase version of the PointFinder results file if it exists
    if [ -f "~{samplename}_PointFinder_results.tsv" ]; then
      awk -F '\t' 'BEGIN{OFS="\t"} { $4=toupper($4) } 1' "~{samplename}_PointFinder_results.tsv" > "~{samplename}_PointFinder_results.uppercase.tsv"
    fi

    # strip off 18 lines from top of file (18th line is the header with the columns: antimicrobial, class, WGS-predicted phenotype, Match, Genetic Background), and convert all letters in first column (antibiotic) to uppercase for readability of output string
    tail +18 ~{samplename}_pheno_table.tsv |  awk -F '\t' 'BEGIN{OFS="\t"} { $1=toupper($1) } 1' > ~{samplename}_pheno_table.headerless.uppercase.tsv
    
    # if column 3 shows 'Resistant', then print list of drugs followed by the genes/point mutations responsible - alphabetized & whitespace trimmed w/ xargs
    awk -F '\t' 'BEGIN{OFS=":"; ORS="; "} { if($3 == "Resistant") {print $1,$5}}' ~{samplename}_pheno_table.headerless.uppercase.tsv | sed 's/..$//' | tr ';' '\n' | sort | tr '\n' ';' | xargs > RESFINDER_PREDICTED_PHENO_RESISTANCE.txt


    # check for XDR Shigella status, based on CDC definition here: https://emergency.cdc.gov/han/2023/han00486.asp
    # requirements: organism input (i.e. gambit_predicted_taxon) must contain "Shigella"
    # predicted resistance to antimicrobials must include ALL: "ampicillin", "azithromycin", "ceftriaxone", "ciprofloxacin", "sulfamethoxazole", and "trimethoprim"
    required_drugs=(ampicillin azithromycin ceftriaxone ciprofloxacin sulfamethoxazole trimethoprim)
    found_all=true

    for drug in "${required_drugs[@]}"; do
      if ! grep -qi "$drug" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
        found_all=false
        break
      fi
    done

    if [[ "~{organism}" == *"Shigella"* ]]; then
      if $found_all; then
        echo "XDR Shigella based on predicted resistance to ampicillin, azithromycin, ceftriaxone, ciprofloxacin, sulfamethoxazole, and trimethoprim. Please verify by reviewing ~{samplename}_pheno_table.tsv and ~{samplename}_ResFinder_results_tab.tsv"
        echo "Potentially XDR Shigella" > RESFINDER_PREDICTED_XDR_SHIGELLA.txt
      else
        echo "Not XDR Shigella" | tee RESFINDER_PREDICTED_XDR_SHIGELLA.txt
      fi
    else
      echo 'Either the input gambit_predicted_taxon does not contain "Shigella" or the user did not supply the organism as an input string to the workflow.'
      echo "Skipping XDR Shigella check."
      echo "Not Shigella based on gambit_predicted_taxon or user input" | tee RESFINDER_PREDICTED_XDR_SHIGELLA.txt
    fi
        
    # function to set output strings for "Resistance" or "No resistance predicted" for drug
    check_resistance() {
      local drug="$1"
      local outfile="$2"
      local tsv="~{samplename}_pheno_table.headerless.uppercase.tsv"
      if grep -qi "$drug" RESFINDER_PREDICTED_PHENO_RESISTANCE.txt; then
        awk -F '\t' -v d="$drug" 'BEGIN{OFS=":";} { if($1 == toupper(d)) {print "Resistance (" $1,$5 ")"}}' "$tsv" > "$outfile"
      else
        echo "No resistance predicted" > "$outfile"
      fi
    }

    # List of drugs and output files
    declare -A drug_to_file=(
      ["ampicillin"]="RESFINDER_PREDICTED_RESISTANCE_AMP.txt"
      ["azithromycin"]="RESFINDER_PREDICTED_RESISTANCE_AZM.txt"
      ["ceftriaxone"]="RESFINDER_PREDICTED_RESISTANCE_AXO.txt"
      ["ciprofloxacin"]="RESFINDER_PREDICTED_RESISTANCE_CIP.txt"
      ["sulfamethoxazole"]="RESFINDER_PREDICTED_RESISTANCE_SMX.txt"
      ["trimethoprim"]="RESFINDER_PREDICTED_RESISTANCE_TMP.txt"
    )

    # Loop through drugs
    for drug in "${!drug_to_file[@]}"; do
      check_resistance "$drug" "${drug_to_file[$drug]}"
    done

    # quinolone resistance detection
    declare -A drug_to_variable=(
      ["CIPROFLOXACIN"]="ciprofloxacin"
      ["FLUOROQUINOLONE"]="fluoroquinolone"
      ["NALIDIXIC ACID"]="nalidixic"
      ["UNKNOWN QUINOLONE"]="unknown_quinolone"
    )
    declare -A resfinder_results
    declare -A pointfinder_results
    declare -A combined_results

    # extract drug information from the pheno table
    for drug in "${!drug_to_variable[@]}"; do
      variable="${drug_to_variable[$drug]}"
      # extract resfinder results
      resfinder_results["$variable"]=$(grep "$drug" ~{samplename}_pheno_table.headerless.uppercase.tsv | awk -F '\t' 'BEGIN{OFS="\t"} { if($3 == "Resistant") {print $5}}' | sed 's/, /\n/g')
      echo "${resfinder_results["$variable"]}" >> quinolone_resistance_candidates.tsv

      # grab pointfinder results if available
      if [ -f "~{samplename}_PointFinder_results.uppercase.tsv" ]; then
        pointfinder_results["$variable"]=$(grep -iE "$drug" ~{samplename}_PointFinder_results.uppercase.tsv | awk -F '\t' 'BEGIN{OFS=""} { split($1,a," "); print a[1] " (" a[2] ")" }')
        echo "${pointfinder_results["$variable"]}" >> quinolone_resistance_candidates.tsv
      fi

      # combine resfinder and pointfinder results
      combined_results["$variable"]=$(printf "%s\n%s\n" "${resfinder_results["$variable"]}" "${pointfinder_results["$variable"]}" | sort -u | grep -v '^$' | paste -d, -s | sed 's/,/, /g')

    done

    # make output string -- prefix with "Resistance" if there are any results, if there are none say "No resistance predicted" instead
    {
      for drug in "${!drug_to_variable[@]}"; do
        variable="${drug_to_variable[$drug]}"
        if [[ -n "${combined_results[$variable]}" ]]; then
          echo "$drug:${combined_results[$variable]}"
        fi
      done
    } | sort | paste -d';' -s > RESFINDER_PREDICTED_RESISTANCE_Q.txt

    # If the output file is empty, write "No resistance predicted"; otherwise prefix and suffix with Resistance (...)
    if [[ ! -s RESFINDER_PREDICTED_RESISTANCE_Q.txt ]]; then
      echo "No resistance predicted" > RESFINDER_PREDICTED_RESISTANCE_Q.txt
    else 
      # add prefix and suffix
      sed -i '1s/^/Resistance (/;1s/$/)/' RESFINDER_PREDICTED_RESISTANCE_Q.txt
    fi
    # add up the number of mechanisms
    sort -u quinolone_resistance_candidates.tsv | grep -v '^$' | wc -l > RESFINDER_PREDICTED_RESISTANCE_Q_COUNT.txt

  >>>
  output {
    File resfinder_pheno_table = "~{samplename}_pheno_table.tsv"
    
    # only if resfinder_organism is set
    File? resfinder_pheno_table_species = "~{samplename}_pheno_table_species.tsv"

    # these files need to be optional in the case where acquired is false
    File? resfinder_hit_in_genome_seq = "~{samplename}_ResFinder_Hit_in_genome_seq.fsa"
    File? resfinder_resistance_gene_seq = "~{samplename}_ResFinder_Resistance_gene_seq.fsa"
    File? resfinder_results_tab = "~{samplename}_ResFinder_results_tab.tsv"

    # only appear if pointfinder = true
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

    String resfinder_predicted_resistance_quinolone = read_string("RESFINDER_PREDICTED_RESISTANCE_Q.txt")
    Int resfinder_predicted_resistance_quinolone_mechanisms = read_string("RESFINDER_PREDICTED_RESISTANCE_Q_COUNT.txt")
    
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
