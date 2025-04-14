version 1.0

task metaspades_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/spades:4.1.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? kmers
    String? metaspades_opts
    Int phred_offset = 33
  }
  command <<<
    # fail hard
    set -euo pipefail

    # get version
    spades.py --meta --version | sed -Ee "s/SPAdes genome assembler ([^ ]+).*/\1/" | tee VERSION

    spades.py \
      --meta \
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

task metaviralspades_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/spades:4.1.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? kmers
    String? metaviralspades_opts
    Int phred_offset = 33
  }
  command <<<
    # fail hard
    set -euo pipefail

    # get version
    spades.py --metaviral --version | sed -Ee "s/SPAdes genome assembler ([^ ]+).*/\1/" | tee VERSION

    # run metaviralspades
    spades.py \
      --metaviral \
      -1 ~{read1_cleaned} \
      -2 ~{read2_cleaned} \
      ~{'-k ' + kmers} \
      -m ~{memory} \
      -t ~{cpu} \
      -o metaspades \
      --phred-offset ~{phred_offset} \
      ~{metaviralspades_opts}

    mv metaspades/contigs.fasta ~{samplename}_contigs.fasta

  >>>
  output {
    File assembly_fasta = "~{samplename}_contigs.fasta"
    String metaviralspades_version = read_string("VERSION")
    String metaviralspades_docker = '~{docker}'
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