version 1.0

task get_fasta_genome_size {
  meta {
    description: "Parse FASTA file to know the length of the genome"
  }
  input {
    File fasta
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/seqkit:2.4.0--h9ee0642_0"
    Int disk_size = 10
    Int cpu = 1
    Int memory = 2
  }
  command <<<
    # from a fasta file, get the length of the genome
    seqkit stats --tabular ~{fasta} > SEQKIT_STDOUT

    if [ ! -s SEQKIT_STDOUT ]; then
      echo "ERROR: seqkit stats failed to generate output!"
      exit 1
    fi

    # verify that STDOUT file only has two lines
    if [ $(cat SEQKIT_STDOUT | wc -l) -ne 2 ]; then
      echo "ERROR: seqkit stats output is not as expected!"
      exit 1
    fi

    # extract the length of the genome - assumes the value is in the fifth column in the last line
    # header: file    format  type    num_seqs        sum_len min_len avg_len max_len
    cat SEQKIT_STDOUT | tail -1 | cut -f 5  | tee GENOME_LENGTH
    
  >>>
  output {
    Int fasta_length = read_int("GENOME_LENGTH")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}