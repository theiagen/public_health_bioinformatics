version 1.0

task ts_mlst {
  meta {
    description: "Torsten Seeman's (TS) automatic MLST calling from assembled contigs"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/mlst:2.23.0-2024-12-31"
    Int disk_size = 50
    Int cpu = 1
    Int memory = 2
    # Parameters
    # --nopath          Strip filename paths from FILE column (default OFF)
    # --scheme [X]      Don't autodetect, force this scheme on all inputs (default '')
    # --minid [n.n]     DNA %identity of full allelle to consider 'similar' [~] (default '95')
    # --mincov [n.n]    DNA %cov to report partial allele at all [?] (default '10')
    # --minscore [n.n]  Minumum score out of 100 to match a scheme (when auto --scheme) (default '50')
    Boolean nopath = true
    Boolean run_secondary_scheme = true # If true, will run secondary scheme if primary scheme is ecoli or abaumannii.
    String? scheme
    String? taxonomy
    Float? min_percent_identity
    Float? min_percent_coverage
    Float? minscore
  }
  command <<< 
    set -euo pipefail

    echo $(mlst --version 2>&1) | sed 's/mlst //' | tee VERSION
    
    #create output header
    echo -e "Filename\tPubMLST_Scheme_name\tSequence_Type_(ST)\tAllele_IDs" > ~{samplename}_ts_mlst.tsv

    # If taxon is E. coli, common mis-characterizations will be excluded from the scheme list
    if [[ "~{taxonomy}" == "Escherichia" || "~{taxonomy}" == "Escherichia coli" || "~{taxonomy}" = "Escherichia_coli" ]]; then
      echo "Taxonomy is ~{taxonomy}, excluding common mis-characterizations for E. coli from the scheme list; aeromonas, cfreundii, senterica"
      mlst \
        --threads ~{cpu} \
        ~{true="--nopath" false="" nopath} \
        ~{'--scheme ' + scheme} \
        ~{'--minid ' + min_percent_identity} \
        ~{'--mincov ' + min_percent_coverage} \
        ~{'--minscore ' + minscore} \
        --exclude 'aeromonas,cfreundii,senterica'
        --novel ~{samplename}_novel_mlst_alleles.fasta \
        ~{assembly} \
        >> ~{samplename}_ts_mlst.tsv
    else 
      mlst \
        --threads ~{cpu} \
        ~{true="--nopath" false="" nopath} \
        ~{'--scheme ' + scheme} \
        ~{'--minid ' + min_percent_identity} \
        ~{'--mincov ' + min_percent_coverage} \
        ~{'--minscore ' + minscore} \
        --novel ~{samplename}_novel_mlst_alleles.fasta \
        ~{assembly} \
        >> ~{samplename}_ts_mlst.tsv
    fi

    if [[ ~{run_secondary_scheme} == true ]]; then
      #create output header
      echo -e "Filename\tPubMLST_Scheme_name\tSequence_Type_(ST)\tAllele_IDs" > ~{samplename}_ts_mlst_secondary_scheme.tsv
      mlst --list > SCHEME_LIST
      scheme=$(head -n 2 ~{samplename}_ts_mlst.tsv | tail -n 1 | cut -f2)
      echo "Scheme Initial Run: $scheme"
      secondary_scheme=$(if [[ "$scheme" == "ecoli_achtman_4" ]]; then
          echo "ecoli"
        elif [[ "$scheme" == "ecoli" ]]; then
          echo "ecoli_achtman_4"
        elif [[ "$scheme" == "abaumannii" ]]; then
          echo "abaumannii_2"
        elif [[ "$scheme" == "abaumannii_2" ]]; then
          echo "abaumannii"
        else
          echo "na"
        fi
      )
      echo "Secondary Scheme: $secondary_scheme"
      if grep -q "$secondary_scheme" SCHEME_LIST; then
        mlst \
          --threads ~{cpu} \
          ~{true="--nopath" false="" nopath} \
          --scheme ${secondary_scheme} \
          ~{'--minid ' + min_percent_identity} \
          ~{'--mincov ' + min_percent_coverage} \
          ~{'--minscore ' + minscore} \
          --novel ~{samplename}_novel_mlst_alleles_secondary_scheme.fasta \
          ~{assembly} \
          >> ~{samplename}_ts_mlst_secondary_scheme.tsv
        
        # parse ts mlst tsv for relevant outputs
        # if output TSV only contains one line (header line); no ST predicted
        if [ $(wc -l ~{samplename}_ts_mlst_secondary_scheme.tsv | awk '{ print $1 }') -eq 1 ]; then
          predicted_mlst_secondary="No ST predicted"
          pubmlst_scheme_secondary="NA"
        # else, TSV has more than one line, so parse outputs
        else
          pubmlst_scheme_secondary="$(cut -f2 ~{samplename}_ts_mlst_secondary_scheme.tsv | tail -n 1)"
          predicted_mlst_secondary="ST$(cut -f3 ~{samplename}_ts_mlst_secondary_scheme.tsv | tail -n 1)"
          # allelic_profile: take second line of output TSV; cut to take 4th column and beyond; replace tabs with commas
          allelic_profile_secondary="$(cut -f 4- ~{samplename}_ts_mlst_secondary_scheme.tsv | tail -n 1 | sed -e 's|\t|,|g')"
          if [ "$pubmlst_scheme_secondary" == "-" ]; then
            predicted_mlst_secondary="No ST predicted"
            pubmlst_scheme_secondary="NA"
          else
            if [ "$predicted_mlst_secondary" == "ST-" ]; then
              predicted_mlst_secondary="No ST predicted"
            fi
          fi
        fi

        echo "$predicted_mlst_secondary" | tee PREDICTED_SECONDARY_MLST
        echo "$pubmlst_scheme_secondary" | tee PUBMLST_SECONDARY_SCHEME
        echo "$allelic_profile_secondary" | tee SECONDARY_ALLELIC_PROFILE.txt
      else
        echo "Secondary scheme $secondary_scheme not found in scheme list: $scheme_list"
        echo "NA" | tee PREDICTED_SECONDARY_MLST
        echo "NA" | tee PUBMLST_SECONDARY_SCHEME
        echo "NA" | tee SECONDARY_ALLELIC_PROFILE.txt
      fi 
    else
      echo "Secondary scheme not run, as run_secondary_scheme is false."
      echo "NA" | tee PREDICTED_SECONDARY_MLST
      echo "NA" | tee PUBMLST_SECONDARY_SCHEME
      echo "NA" | tee SECONDARY_ALLELIC_PROFILE.txt
    fi

    # Will run on primary scheme results
    # parse ts mlst tsv for relevant outputs
    # if output TSV only contains one line (header line); no ST predicted
    if [ $(wc -l ~{samplename}_ts_mlst.tsv | awk '{ print $1 }') -eq 1 ]; then
      predicted_mlst="No ST predicted"
      pubmlst_scheme="NA"
    # else, TSV has more than one line, so parse outputs
    else
      pubmlst_scheme="$(cut -f2 ~{samplename}_ts_mlst.tsv | tail -n 1)"
      predicted_mlst="ST$(cut -f3 ~{samplename}_ts_mlst.tsv | tail -n 1)"
      # allelic_profile: take second line of output TSV; cut to take 4th column and beyond; replace tabs with commas
      allelic_profile="$(cut -f 4- ~{samplename}_ts_mlst.tsv | tail -n 1 | sed -e 's|\t|,|g')"
      if [ "$pubmlst_scheme" == "-" ]; then
        predicted_mlst="No ST predicted"
        pubmlst_scheme="NA"
      else
        if [ "$predicted_mlst" == "ST-" ]; then
          predicted_mlst="No ST predicted"
        fi
      fi  
    fi
        
    echo "$predicted_mlst" | tee PREDICTED_MLST
    echo "$pubmlst_scheme" | tee PUBMLST_SCHEME
    echo "$allelic_profile" | tee ALLELIC_PROFILE.txt

    # Concatenates the secondary scheme results if it exists
    if [[ -f "~{samplename}_ts_mlst_secondary_scheme.tsv" ]]; then
      awk 'NR>1' ~{samplename}_ts_mlst_secondary_scheme.tsv > ~{samplename}_ts_mlst_secondary_scheme_no_header.tsv
      cat ~{samplename}_ts_mlst_secondary_scheme_no_header.tsv >> ~{samplename}_ts_mlst.tsv
    fi
  >>>
  output {
    File ts_mlst_results = "~{samplename}_ts_mlst.tsv"
    String ts_mlst_predicted_st = read_string("PREDICTED_MLST")
    String ts_mlst_pubmlst_scheme = read_string("PUBMLST_SCHEME")
    String ts_mlst_allelic_profile = read_string("ALLELIC_PROFILE.txt")
    File? ts_mlst_novel_alleles = "~{samplename}_novel_mlst_alleles.fasta"
    # Only present if secondary scheme was run
    String? ts_mlst_predicted_secondary_st = read_string("PREDICTED_SECONDARY_MLST")
    String? ts_mlst_pubmlst_secondary_scheme = read_string("PUBMLST_SECONDARY_SCHEME") 
    String? ts_mlst_secondary_allelic_profile = read_string("SECONDARY_ALLELIC_PROFILE.txt")
    File? ts_mlst_secondary_novel_alleles = "~{samplename}_novel_mlst_alleles_secondary_scheme.fasta"
    String ts_mlst_version = read_string("VERSION")
    String ts_mlst_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}
