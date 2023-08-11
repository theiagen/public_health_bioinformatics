version 1.0

task resfinder {
  input {
    File assembly # Input fasta file
    String samplename
    String? organism # Species in the sample, species should be entered with their full scientific names (e.g. "escherichia coli"), using quotation marks
    Boolean acquired = true # Run resfinder for acquired resistance genes
    Float min_cov = 0.6 # Minimum (breadth-of) coverage of ResFinder
    Float min_id = 0.9 # Threshold for identity of ResFinder
    Boolean call_PointFinder = false # Run pointfinder for chromosomal mutations
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
        ~{true="--point" false="" call_PointFinder} 
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

  >>>
  output {
    File resfinder_pheno_table = "~{samplename}_pheno_table.tsv"
    File? resfinder_pheno_table_species = "~{samplename}_pheno_table_species.tsv"
    File resfinder_hit_in_genome_seq = "~{samplename}_ResFinder_Hit_in_genome_seq.fsa"
    File resfinder_resistance_gene_seq = "~{samplename}_ResFinder_Resistance_gene_seq.fsa"
    File resfinder_results_tab = "~{samplename}_ResFinder_results_tab.tsv"
    File? pointfinder_pheno_table = "~{samplename}_PointFinder_prediction.tsv"
    File? pointfinder_results = "~{samplename}_PointFinder_results.tsv"
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
