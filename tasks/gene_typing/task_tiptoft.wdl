version 1.0

task tiptoft {
  input {
    String samplename
    File read1 # intended for ONT data only
    Int disk_size = 100
    Int cpu = 2
    Int? kmer_size # default is 13
    Int? max_gap # default is 3
    Int? margin # default is 10
    Int? min_block_size # default is 130
    Int? min_fasta_hits # default is 10
    Int? min_perc_coverage # default is 85
    Int? min_kmers_for_onex_pass # default is 10
    String docker = "staphb/tiptoft:1.0.2"
  }
  command <<<
    # capture version
    tiptoft --version | tee VERSION

    # run TipToft on FASTQ file; leave all options as default, but allow user to adjust if necessary
    tiptoft \
      ~{read1} \
      ~{'--kmer ' + kmer_size} \
      ~{'--max_gap ' + max_gap} \
      ~{'--margin ' + margin} \
      ~{'--min_block_size ' + min_block_size} \
      ~{'--min_fasta_hits ' + min_fasta_hits} \
      ~{'--min_perc_coverage ' + min_perc_coverage} \
      ~{'--min_kmers_for_onex_pass ' + min_kmers_for_onex_pass} \
      --filtered_reads_file ~{samplename}.tiptoft.fastq \
      --output_file ~{samplename}.tiptoft.tsv

    # gzip output filtered FASTQ; output filename will be ~{samplename}.tiptoft.fastq.gz
    # could change this to pigz if we include it in the docker image
    gzip ~{samplename}.tiptoft.fastq

    # parse tiptoft output TSV
    # parse out gene names into list of strings, comma-separated, final comma at end removed by sed
    plasmid_replicon_genes=$(awk -F '\t' '{ print $1 }' ~{samplename}.tiptoft.tsv | tail -n+2 | tr '\n' ',' | sed 's/.$//')
    echo "${plasmid_replicon_genes}" | tee PLASMID_REPLICON_GENES.txt

  >>>
  output {
    File tiptoft_tsv = "~{samplename}.tiptoft.tsv"
    File tiptoft_plasmid_replicon_fastq = "~{samplename}.tiptoft.fastq.gz"
    String plasmid_replicon_genes = read_string("PLASMID_REPLICON_GENES.txt")
    String tiptoft_version = read_string("VERSION")
    String tiptoft_docker = docker
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}
