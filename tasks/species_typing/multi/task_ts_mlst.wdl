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
    String? scheme
    String? taxonomy
    Float? min_percent_identity
    Float? min_percent_coverage
    Float? minscore
  }
  command <<<
    echo $(mlst --version 2>&1) | sed 's/mlst //' | tee VERSION
    
    #create output header
    echo -e "Filename\tPubMLST_Scheme_name\tSequence_Type_(ST)\tAllele_IDs" > ~{samplename}_ts_mlst.tsv

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

    scheme=$(head -n 2 ~{samplename}_ts_mlst.tsv | tail -n 1 | cut -f2)
    echo "Scheme Initial Run: $scheme"

    if [[ "$scheme" == "ecoli" || "$scheme" == "ecoli_2" || "$scheme" == "abaumannii" || "$scheme" == "abaumannii_2" ]]; then
      mv ~{samplename}_ts_mlst.tsv ~{samplename}_1.tsv
      secondary_scheme=$(if [[ "$scheme" == *_2 ]]; then echo "${scheme%_2}"; else echo "${scheme}_2"; fi)

      mlst \
      --threads ~{cpu} \
      ~{true="--nopath" false="" nopath} \
      --scheme $secondary_scheme \
      ~{'--minid ' + min_percent_identity} \
      ~{'--mincov ' + min_percent_coverage} \
      ~{'--minscore ' + minscore} \
      --novel ~{samplename}_novel_mlst_alleles_${secondary_scheme}.fasta \
      ~{assembly} \
      >> ~{samplename}_2.tsv

      #create output header

      cat ~{samplename}_1.tsv ~{samplename}_2.tsv >> ~{samplename}_ts_mlst.tsv
      cat ~{samplename}_novel_mlst_alleles_${secondary_scheme}.fasta ~{samplename}_novel_mlst_alleles.fasta > ~{samplename}_novel_mlst_alleles.fasta

    elif [[ "~{taxonomy}" == "Escherichia" || "~{taxonomy}" == "Escherichia coli" ]]; then
      if [[ "$scheme" == "aeromonas" || "$scheme" == "cfreundii" ||  "$scheme" == "senterica" ]]; then
          echo "Taxonomy is reported as Escherichia, but scheme is not. Running ecoli schemes."

          mlst \
          --threads ~{cpu} \
          ~{true="--nopath" false="" nopath} \
          --scheme ecoli \
          ~{'--minid ' + min_percent_identity} \
          ~{'--mincov ' + min_percent_coverage} \
          ~{'--minscore ' + minscore} \
          --novel ~{samplename}_novel_mlst_alleles_ecoli.fasta \
          ~{assembly} \
          >> ~{samplename}_1.tsv

          mlst \
          --threads ~{cpu} \
          ~{true="--nopath" false="" nopath} \
          --scheme ecoli_2 \
          ~{'--minid ' + min_percent_identity} \
          ~{'--mincov ' + min_percent_coverage} \
          ~{'--minscore ' + minscore} \
          --novel ~{samplename}_novel_mlst_alleles_ecoli_2.fasta \
          ~{assembly} \
          >> ~{samplename}_2.tsv

          #create output header
          echo -e "Filename\tPubMLST_Scheme_name\tSequence_Type_(ST)\tAllele_IDs" > ~{samplename}_ts_mlst.tsv

          cat ~{samplename}_1.tsv ~{samplename}_2.tsv >> ~{samplename}_ts_mlst.tsv
          cat ~{samplename}_novel_mlst_alleles_ecoli.fasta ~{samplename}_novel_mlst_alleles_ecoli_2.fasta > ~{samplename}_novel_mlst_alleles.fasta
      fi
    fi
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
  >>>
  output {
    File? ts_mlst_results = "~{samplename}_ts_mlst.tsv"
    String ts_mlst_predicted_st = read_string("PREDICTED_MLST")
    String ts_mlst_pubmlst_scheme = read_string("PUBMLST_SCHEME")
    String ts_mlst_allelic_profile = read_string("ALLELIC_PROFILE.txt")
    File? ts_mlst_novel_alleles = "~{samplename}_novel_mlst_alleles.fasta"
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
