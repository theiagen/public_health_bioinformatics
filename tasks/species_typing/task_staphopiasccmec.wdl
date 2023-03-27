version 1.0

task staphopiasccmec {
  meta {
    description: "Primer based SCCmec typing of Staphylococcus aureus genomes"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0"
    Int disk_size = 100
    Int cpu = 1
  }
  command <<<
    # get version
    staphopia-sccmec --version 2>&1 | sed 's/^.*staphopia-sccmec //' | tee VERSION

    # run staphopia-sccmec on input assembly; hamming option OFF; outputs are true/false
    staphopia-sccmec \
      --assembly ~{assembly} > ~{samplename}.staphopia-sccmec.summary.tsv

    # run staphopia-sccmec on input assembly; hamming option ON; outputs are the hamming distance; 0 is exact match
    staphopia-sccmec \
      --hamming \
      --assembly ~{assembly} > ~{samplename}.staphopia-sccmec.hamming.tsv

    # please excuse this ugly bash code below :)

    # parse output summary TSV for true matches
    # look for columns that contain the word "True" and print the column numbers in a list to a file col_headers.txt
    awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "True") print i } }' ~{samplename}.staphopia-sccmec.summary.tsv | tee col_headers.txt

    # use column number list to print column headers (example: IV, mecA, etc.) to a file type.txt
    cat col_headers.txt | while read -r COL_NUMBER; do \
      cut -f "$COL_NUMBER" ~{samplename}.staphopia-sccmec.summary.tsv | head -n 1 >>type.txt
      echo "," >>type.txt
    done

    # remove newlines, remove trailing comma; generate output string of comma separated values
    cat type.txt | tr -d '\n' | sed 's|.$||g' | tee TYPES_AND_MECA.txt
  >>>
  output {
    File staphopiasccmec_results_tsv = "~{samplename}.staphopia-sccmec.summary.tsv"
    File staphopiasccmec_hamming_distance_tsv = "~{samplename}.staphopia-sccmec.hamming.tsv"
    String staphopiasccmec_types_and_mecA_presence = read_string("TYPES_AND_MECA.txt")
    String staphopiasccmec_version = read_string("VERSION")
    String staphopiasccmec_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
