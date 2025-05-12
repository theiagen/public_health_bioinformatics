version 1.0

task spades {
  input {
    File read1
    File? read2
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
    spades.py --version | sed -Ee "s/SPAdes genome assembler ([^ ]+).*/\1/" | tee VERSION

    echo "DEBUG: Running SPAdes"

    if [ -n "~{read2}" ]; then
      spades.py \
        ~{'--' + spades_type} \
        -1 ~{read1} \
        -2 ~{read2} \
        ~{'-k ' + kmers} \
        -m ~{memory} \
        -t ~{cpu} \
        -o spades \
        --phred-offset ~{phred_offset} \
        ~{spades_opts}
    else
      spades.py \
        ~{'--' + spades_type} \
        -s ~{read1} \
        ~{'-k ' + kmers} \
        -m ~{memory} \
        -t ~{cpu} \
        -o spades \
        --phred-offset ~{phred_offset} \
        ~{spades_opts}
    fi

    if [ ~{spades_type} == "metaviral" ]; then
      # if metaviralspades fails, or fails to output a contigs.fasta, we want to report that for falling back
      if [ -f spades/contigs.fasta ]; then
        echo "DEBUG: Metaviralspades successfully identified a complete virus"
        mv spades/contigs.fasta ~{samplename}~{'_' + spades_type + 'spades'}_contigs.fasta
        # Move the contigs.gfa file to the top level directory
        mv spades/assembly_graph_with_scaffolds.gfa ~{samplename}~{'_' + spades_type + 'spades'}_contigs.gfa
        echo "PASS" > STATUS
      else
        echo "DEBUG: Metaviralspades could not identify a complete virus"
        echo "FAIL" > STATUS
      fi
    else
      # all other spades types should fail via the pipefail, so they pass by default
      mv spades/contigs.fasta ~{samplename}~{'_' + spades_type + 'spades'}_contigs.fasta
      # Move the contigs.gfa file to the top level directory
      mv spades/assembly_graph_with_scaffolds.gfa ~{samplename}~{'_' + spades_type + 'spades'}_contigs.gfa
      echo "PASS" > STATUS
    fi
  >>>
  output {
    File assembly_fasta = "~{samplename}~{'_' + spades_type + 'spades'}_contigs.fasta"
    File? assembly_gfa = "~{samplename}~{'_' + spades_type + 'spades'}_contigs.gfa"
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