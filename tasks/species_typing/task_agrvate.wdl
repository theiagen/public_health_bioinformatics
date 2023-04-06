version 1.0

task agrvate {
  meta {
    description: "Rapid identification of Staphylococcus aureus agr locus type and agr operon variants."
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/agrvate:1.0.2--hdfd78af_0"
    Int disk_size = 50
    Int cpu = 1

    # Parameters
    # --typing_only    agr typing only. Skips agr operon extraction and frameshift detection
    Boolean typing_only = false
  }
  command <<<
    # get version info
    agrvate -v 2>&1 | sed 's/agrvate v//;' | tee VERSION

    # run agrvate on assembly; usearch not available in biocontainer, cannot use that option
    # using -m flag for mummer frameshift detection since usearch is not available
    agrvate \
        ~{true="--typing-only" false="" typing_only} \
        -i ~{assembly} \
        -m 

    # agrvate names output directory and file based on name of .fasta file, so <prefix>.fasta as input results in <prefix>-results/ outdir
    # and results in <prefix>-results/<prefix>-summary.tab files 
    basename=$(basename ~{assembly})
    # strip off anything after the period
    fasta_prefix=${basename%.*}

    # rename outputs summary TSV to include samplename
    mv -v "${fasta_prefix}-results/${fasta_prefix}-summary.tab" ~{samplename}.agrvate.tsv

    # parse output summary TSV
    cut -f 2 ~{samplename}.agrvate.tsv | tail -n 1 | tee AGR_GROUP
    cut -f 3 ~{samplename}.agrvate.tsv | tail -n 1 | tee AGR_MATCH_SCORE
    cut -f 4 ~{samplename}.agrvate.tsv | tail -n 1 | tee AGR_CANONICAL
    cut -f 5 ~{samplename}.agrvate.tsv | tail -n 1 | tee AGR_MULTIPLE
    cut -f 6 ~{samplename}.agrvate.tsv | tail -n 1 | tee AGR_NUM_FRAMESHIFTS

    # edit output string AGR_CANONICAL to be more informative: https://github.com/VishnuRaghuram94/AgrVATE#results
    if [[ $(cat AGR_CANONICAL) == 1 ]]; then
      echo "1. canonical agrD" >AGR_CANONICAL
    elif [[ $(cat AGR_CANONICAL) == 0 ]]; then
      echo "0. non-canonical agrD" >AGR_CANONICAL
    elif [[ $(cat AGR_CANONICAL) == "u" ]]; then
      echo "u. unknown agrD" >AGR_CANNONICAL
    else 
      echo "result unrecognized, please see summary agrvate TSV file" >AGR_CANONICAL
    fi

    # edit output string AGR_MULTIPLE to be more informative: https://github.com/VishnuRaghuram94/AgrVATE#results
    if [[ $(cat AGR_MULTIPLE) == "s" ]]; then
      echo "s. single agr group found" >AGR_MULTIPLE
    elif [[ $(cat AGR_MULTIPLE) == "m" ]]; then
      echo "m. multiple agr groups found" >AGR_MULTIPLE
    elif [[ $(cat AGR_MULTIPLE) == "u" ]]; then
      echo "u. unknown agr groups found" >AGR_MULTIPLE
    else 
      echo "result unrecognized, please see summary agrvate TSV file" >AGR_MULTIPLE
    fi

    # if AGR_NUM_FRAMESHIFTS is unknown, edit output string AGR_NUM_FRAMESHIFTS to be more informative, otherwise keep set to a number: https://github.com/VishnuRaghuram94/AgrVATE#results
    if [[ $(cat AGR_NUM_FRAMESHIFTS) == "u" ]]; then
      echo "u or unknown; agr operon not extracted" >AGR_NUM_FRAMESHIFTS
    fi

    # create tarball of all output files
    tar -czvf ~{samplename}.agrvate.tar.gz "${fasta_prefix}-results/"
  >>>
  output {
    File agrvate_summary = "~{samplename}.agrvate.tsv"
    File agrvate_results = "~{samplename}.agrvate.tar.gz"
    String agrvate_agr_group = read_string("AGR_GROUP")
    String agrvate_agr_match_score = read_string("AGR_MATCH_SCORE")
    String agrvate_agr_canonical = read_string("AGR_CANONICAL")
    String agrvate_agr_multiple = read_string("AGR_MULTIPLE")
    String agrvate_agr_num_frameshifts = read_string("AGR_NUM_FRAMESHIFTS")
    String agrvate_version = read_string("VERSION")
    String agrvate_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}
