version 1.0

task metaspades_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/spades:3.12.0--h9ee0642_3"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? kmers
    String? metaspades_opts
    Int phred_offset = 33
  }
  command <<<
    metaspades.py --version | head -1 | cut -d ' ' -f 2 | tee VERSION
    metaspades.py \
      -1 ~{read1_cleaned} \
      -2 ~{read2_cleaned} \
      ~{'-k ' + kmers} \
      -m ~{memory} \
      -t ~{cpu} \
      -o metaspades \
      --phred-offset ~{phred_offset} \
      ~{metaspades_opts}

    mv metaspades/contigs.fasta ~{samplename}_contigs.fasta

  >>>
  output {
    File assembly_fasta = "~{samplename}_contigs.fasta"
    String metaspades_version = read_string("VERSION")
    String metaspades_docker = '~{docker}'
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

