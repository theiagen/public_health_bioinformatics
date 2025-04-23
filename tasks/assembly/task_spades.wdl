version 1.0

task spades {
  input {
    File read1_cleaned
    File? read2_cleaned
    String samplename
    String? spades_type
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/spades:4.1.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16
    String? kmers
    String? spades_opts
    Int phred_offset = 33
  }
  command <<<
    # fail hard if not metaviralspades
    if [ ! ~{spades_type} == "metaviral" ]; then
      set -euo pipefail
    fi

    # get version
    spades.py ${spades_call} --version | sed -Ee "s/SPAdes genome assembler ([^ ]+).*/\1/" | tee VERSION

    echo "DEBUG: Running SPAdes"
    spades.py \
      ~{'--' + spades_type} \
      -1 ~{read1_cleaned} \
      ~{'-2 ' + read2_cleaned} \
      ~{'-k ' + kmers} \
      -m ~{memory} \
      -t ~{cpu} \
      -o spades \
      --phred-offset ~{phred_offset} \
      ~{spades_opts}

    if [ ~{spades_type} == "metaviral" ]; then
      # if metaviralspades fails, or fails to output a contigs.fasta, we want to report that for falling back
      if [ -f spades/contigs.fasta ]; then
        echo "DEBUG: Metaviralspades successfully identified a complete virus"
        mv spades/contigs.fasta ~{samplename}~{'_' + spades_type}_contigs.fasta
        echo "PASS" > STATUS
      else
        echo "DEBUG: Metaviralspades could not identify a complete virus"
        echo "FAIL" > STATUS
      fi
    else
      # all other spades types should fail via the pipefail, so they pass by default
      mv spades/contigs.fasta ~{samplename}~{'_' + spades_type}_contigs.fasta
      echo "PASS" > STATUS
    fi
  >>>
  output {
    File assembly_fasta = "~{samplename}~{' ' + spades_type}_contigs.fasta"
    String spades_status = read_string("STATUS")
    String spades_version = read_string("VERSION")
    String spades_docker = '~{docker}'
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