version 1.0

task animummer {
  input {
    File assembly
    String samplename
    File? ref_genome
    Float mash_filter = 0.9
    String docker="staphb/mummer:4.0.0-rgdv2"
  }
  command <<<
    # capture and version
    mummer --version | tee MUMMER_VERSION

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
      echo "0.0" > ANI_HIGHEST_PERCENT_BASES_ALIGNED
      echo "0.0" > ANI_HIGHEST_PERCENT
      echo "ANI skipped due to high genetic divergence from reference genomes" > ANI_TOP_SPECIES_MATCH
    # if output TSV has greater than 1 lines, then parse for appropriate outputs
    else
      ## parse out highest percentBases aligned
      cut -f 5 ~{samplename}.ani-mummer.out.tsv | sort -nr | head -n 1 | tee ANI_HIGHEST_PERCENT_BASES_ALIGNED
      echo "highest percent bases aligned is: $(cat ANI_HIGHEST_PERCENT_BASES_ALIGNED)"

      ## parse out ANI value using highest percentBases aligned value
      grep "$(cat ANI_HIGHEST_PERCENT_BASES_ALIGNED)" ~{samplename}.ani-mummer.out.tsv | cut -f 3 | tee ANI_HIGHEST_PERCENT
      echo "ANI value is: $(cat ANI_HIGHEST_PERCENT)"

      # have to separate out results for ani_top_species match because user-defined reference genome FASTAs will not be named as they are in RGDv2
      if [[ -z "~{ref_genome}" ]]; then
        ### ref genome is not user-defined, using RGDv2 and FASTA filenames ###
        # Parse out species name from reference fasta filename
        # use percent bases aligned to pull relevant line, cut down to query and ref fasta filenames, sed to remove your query filename, xargs to remove whitespaces & stuff
        # cut on periods to pull out genus_species (in future this will inlcude lineages for Listeria and other sub-species designations)
        # have to create assembly_file_basename bash variable since output TSV does not include full path to assembly file, only filename
        assembly_file_basename=$(basename ~{assembly})
        grep "$(cat ANI_HIGHEST_PERCENT)" ~{samplename}.ani-mummer.out.tsv | cut -f 1,2 | sed "s|${assembly_file_basename}||g" | xargs | cut -d '.' -f 3 | tee ANI_TOP_SPECIES_MATCH
        echo "ANI top species match is: $(cat ANI_TOP_SPECIES_MATCH)"
      else 
        # User specified a reference genome, use fasta filename as output string
        basename "${REF_GENOME}" > ANI_TOP_SPECIES_MATCH
        echo "Reference genome used for ANI is: ${REF_GENOME}"
      fi
    fi
    
  >>>
  output {
    Float ani_highest_percent = read_float("ANI_HIGHEST_PERCENT")
    Float ani_highest_percent_bases_aligned = read_float("ANI_HIGHEST_PERCENT_BASES_ALIGNED")
    File ani_output_tsv = "~{samplename}.ani-mummer.out.tsv"
    String ani_top_species_match = read_string("ANI_TOP_SPECIES_MATCH")
    String ani_mummer_version = read_string("MUMMER_VERSION")
  }
  runtime {
    docker:  "~{docker}"
    memory:  "8 GB"
    cpu:   4
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}