version 1.0

task animummer {
  input {
    File assembly
    String samplename
    File? ref_genome
    Float mash_filter = 0.9
    # these 2 thresholds were set as they are used by CDC enterics lab/PulseNet for ANI thresholds
    Float ani_threshold = 80.0
    Float percent_bases_aligned_threshold = 70.0
    String docker= "us-docker.pkg.dev/general-theiagen/staphb/mummer:4.0.0-rgdv2"
    Int cpus = 4
    Int memory = 8
    Int disk_size = 100
  }
  command <<<
    # capture and version
    mummer --version | tee MUMMER_VERSION.txt

    # set the reference genome
    # if not defined by user, then use all 43 genomes in RGDv2
    if [[ -z "~{ref_genome}" ]]; then
      # ref genome is not defined. default to RGDv2
      # BASH variable
      REF_GENOME="$(ls /RGDv2/*.fasta)"
      echo "user did not define a reference genome, defaulting to 43 genomes in RGDv2"
      echo "REF_GENOME is set to: ${REF_GENOME}"
    else 
      echo "User specified a reference genome, will use this instead of RGDv2"
      REF_GENOME="~{ref_genome}"
      echo "REF_GENOME is set to: ${REF_GENOME}"
    fi

    # call Lee's ani-m.pl script and compare query genome against reference genome
    # first does a mash check on relatedness between 2 genomes. If greater than mash_filter, then run dnadiff
    # --symmetric flag runs ANI on query vs. ref; followed by ref vs. query
    ani-m.pl --symmetric \
             --mash-filter ~{mash_filter} \
             ~{assembly} \
             ${REF_GENOME} | tee ~{samplename}.ani-mummer.out.tsv

    # CHECK FOR A NEARLY BLANK TSV (ONLY HEADER LINE), mean sample did not surpass mash-filter and thus no ANI was run
    LINE_COUNT_OUTPUT_TSV=$(wc -l ~{samplename}.ani-mummer.out.tsv | cut -d ' ' -f 1)
    echo "Number of lines in output TSV is: ${LINE_COUNT_OUTPUT_TSV}"
    if [[ ${LINE_COUNT_OUTPUT_TSV} -eq 1 ]]; then
      echo "~{samplename} did not surpass the minimum mash genetic distance filter, thus ANI was not performed"
      echo "The output TSV only contains the header line"
      # set output variables as 0s or descriptive strings
      echo "0.0" > ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt
      echo "0.0" > ANI_HIGHEST_PERCENT.txt
      echo "ANI skipped due to high genetic divergence from reference genomes" > ANI_TOP_SPECIES_MATCH.txt
    # if output TSV has greater than 1 lines, then parse for appropriate outputs
    else
      ## parse out highest percentBases aligned
      cut -f 5 ~{samplename}.ani-mummer.out.tsv | sort -nr | head -n 1 | tee ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt
      echo "highest percent bases aligned is: $(cat ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt)"
      ANI_HIGHEST_PERCENT_BASES_ALIGNED=$(cat ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt)

      ## parse out ANI value using highest percentBases aligned value
      grep "$(cat ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt)" ~{samplename}.ani-mummer.out.tsv | cut -f 3 | tee ANI_HIGHEST_PERCENT.txt
      echo "Highest ANI value is: $(cat ANI_HIGHEST_PERCENT.txt)"
      # set ANI_HIGHEST_PERCENT as a bash variable (float)
      ANI_HIGHEST_PERCENT=$(cat ANI_HIGHEST_PERCENT.txt)


      # have to separate out results for ani_top_species match because user-defined reference genome FASTAs will not be named as they are in RGDv2
      if [[ -z "~{ref_genome}" ]]; then
        ### ref genome is not user-defined, using RGDv2 and FASTA filenames ###
        # Parse out species name from reference fasta filename
        # use percent bases aligned to pull relevant line, cut down to query and ref fasta filenames, sed to remove your query filename, xargs to remove whitespaces & stuff
        # cut on periods to pull out genus_species (in future this will inlcude lineages for Listeria and other sub-species designations)
        # have to create assembly_file_basename bash variable since output TSV does not include full path to assembly file, only filename
        assembly_file_basename=$(basename ~{assembly})
        grep "${ANI_HIGHEST_PERCENT}" ~{samplename}.ani-mummer.out.tsv | cut -f 1,2 | sed "s|${assembly_file_basename}||g" | xargs | cut -d '.' -f 3 | tee ANI_TOP_SPECIES_MATCH.txt
        echo "ANI top species match is: $(cat ANI_TOP_SPECIES_MATCH.txt)"

        # if ANI threshold or percent_bases_aligned_threshold is defined by user (they both are by default), compare to highest ANI value and corresponding percent_bases_aligned value and only output ANI_top_species_match if both thresholds are surpassed
        if [[ -n "~{ani_threshold}" || -n "~{percent_bases_aligned_threshold}" ]]; then
          echo "Comparing user-defined ANI threshold to highest ANI value..."
          # compare ANI_HIGHEST_PERCENT to ani_threshold using awk
          if ! awk "BEGIN{ exit ($ANI_HIGHEST_PERCENT < ~{ani_threshold} )}"; then
            echo "The highest ANI value $ANI_HIGHEST_PERCENT is less than the user-defined ANI threshold of ~{ani_threshold}"
            echo "ANI top species match did not surpass the user-defined ANI threshold of ~{ani_threshold}" > ANI_TOP_SPECIES_MATCH.txt
          # else if: compare percent_bases_aligned_threshold to ANI_HIGHEST_PERCENT_BASES_ALIGNED using awk
          elif ! awk "BEGIN{ exit (${ANI_HIGHEST_PERCENT_BASES_ALIGNED} < ~{percent_bases_aligned_threshold} )}"; then
            echo "The highest ANI percent bases aligned value ${ANI_HIGHEST_PERCENT_BASES_ALIGNED} is less than the user-defined threshold of ~{percent_bases_aligned_threshold}"
            # overwrite ANI_TOP_SPECIES_MATCH.txt when percent_bases_aligned threshold is not surpassed
            echo "ANI top species match did not surpass the user-defined percent bases aligned threshold of ~{percent_bases_aligned_threshold}" > ANI_TOP_SPECIES_MATCH.txt 
          else
            echo "The highest ANI value ${ANI_HIGHEST_PERCENT} is greater than the user-defined threshold ~{ani_threshold}"
            echo "The highest percent bases aligned value ${ANI_HIGHEST_PERCENT_BASES_ALIGNED} is greater than the user-defined threshold ~{percent_bases_aligned_threshold}"
          fi
        fi
      else 
        # User specified a reference genome, use fasta filename as output string
        basename "${REF_GENOME}" > ANI_TOP_SPECIES_MATCH.txt
        echo "Reference genome used for ANI is: ${REF_GENOME}"
      fi
    fi
    
  >>>
  output {
    Float ani_highest_percent = read_float("ANI_HIGHEST_PERCENT.txt")
    Float ani_highest_percent_bases_aligned = read_float("ANI_HIGHEST_PERCENT_BASES_ALIGNED.txt")
    File ani_output_tsv = "~{samplename}.ani-mummer.out.tsv"
    String ani_top_species_match = read_string("ANI_TOP_SPECIES_MATCH.txt")
    String ani_mummer_version = read_string("MUMMER_VERSION.txt")
    String ani_docker = "~{docker}"
  }
  runtime {
    docker:  "~{docker}"
    memory:  "~{memory} GB"
    cpu:  cpus
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}