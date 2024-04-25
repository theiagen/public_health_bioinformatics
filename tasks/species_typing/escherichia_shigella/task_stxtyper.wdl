version 1.0

task stxtyper {
  input {
    File assembly
    String samplename
    String docker = "kapsakcj/stxtyper:8328c4d"
    Int disk_size = 100
    Int cpu = 2
    Int memory = 8
  }
  command <<<
    # capture date
    date | tee DATE

    # capture version info
    stxtyper --version | tee VERSION.txt

    echo "DEBUG: running StxTyper now..."
    # run StxTyper on assembly; may need to add/remove options in the future if they change
    stxtyper \
      --nucleotide ~{assembly} \
      --name ~{samplename} \
      --output ~{samplename}_stxtyper.tsv \
      --log ~{samplename}_stxtyper.log

    # parse output TSV
    echo "DEBUG: Parsing StxTyper output TSV..."

    # check for output file with only 1 line (meaning no hits found); exit cleanly if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -eq 1 ]; then
      echo "No hits found by StxTyper" > stxtyper_hits.txt
      echo "0" > stxtyper_num_hits.txt
      echo "DEBUG: No hits found in StxTyper output TSV. Exiting task with exit code 0 now."
      exit 0
    fi
    
    # check for output file with more than 1 line (meaning hits found); count lines & parse output TSV if so
    if [ "$(wc -l < ~{samplename}_stxtyper.tsv)" -gt 1 ]; then
      echo "Hits found by StxTyper. Counting lines & parsing output TSV now..."
      # count number of lines in output TSV (excluding header)
      wc -l < ~{samplename}_stxtyper.tsv | awk '{print $1-1}' > stxtyper_num_hits.txt
      # remove header line
      sed '1d' ~{samplename}_stxtyper.tsv > ~{samplename}_stxtyper_noheader.tsv
      # parse output TSV and create comma separated list of: (stx_type|operon|target_contig)
      awk -F'\t' '{print "(" $3 "|" $4 "|" $2 ")"}' ~{samplename}_stxtyper_noheader.tsv | paste -sd, - | tee stxtyper_hits.txt
    fi

    echo "DEBUG: Finished parsing StxTyper output TSV."
  >>>
  output {
    File stxtyper_report = "~{samplename}_stxtyper.tsv"
    File stxtyper_log = "~{samplename}_stxtyper.log"
    String stxtyper_docker = docker
    String stxtyper_version = read_string("VERSION.txt")
    String stxtyper_hits = read_string("stxtyper_hits.txt")
    Int stxtyper_num_hits = read_int("stxtyper_num_hits.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
