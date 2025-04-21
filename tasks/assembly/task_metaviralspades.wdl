version 1.0

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

    if [ -f metaspades/contigs.fasta ]; then
      echo "Metaviralspades successfully identified a complete virus"
      mv metaspades/contigs.fasta ~{samplename}_metaviralspades.fasta
      echo "PASS" > STATUS
    else
      echo "Metaviralspades could not idnentify a complete virus"
      echo "FAIL" > STATUS
    fi

  >>>
  output {
    File? assembly_fasta = "~{samplename}_metaviralspades.fasta"
    String metaviralspades_status = read_string("STATUS")
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