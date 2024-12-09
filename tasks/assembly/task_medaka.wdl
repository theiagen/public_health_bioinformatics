version 1.0

task medaka_consensus {
    input {
      File unpolished_fasta
      String samplename
      File read1
      String medaka_model = "r1041_e82_400bps_sup_v5.0.0"
      Int polish_rounds = 1
      Int cpu = 4
      Int memory = 16
      Int disk_size = 100
      String docker = "us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1"
    }
    command <<< 
      set -euo pipefail
      medaka --version | tee MEDAKA_VERSION

      # Initialize the input for polishing with the provided assembly FASTA
        cp ~{unpolished_fasta} polished_input.fasta

        # Perform Medaka polishing for the specified number of rounds
        for i in $(seq 1 ~{polish_rounds}); do
            echo "Starting Medaka polishing round $i"
            polish_dir="polish_round_$i"
            mkdir -p "$polish_dir"

            medaka_consensus \
                -i ~{read1} \
                -d polished_input.fasta \
                -o "$polish_dir" \
                -m ~{medaka_model} \
                -t ~{cpu}

            # Updates the input for the next round
            cp "$polish_dir/consensus.fasta" polished_input.fasta
        done

        # Rename final polished outputs with sample name
        mv polished_input.fasta ~{samplename}.polished.fasta
    >>>
    output {
      File medaka_fasta = "~{samplename}.polished.fasta"
      String medaka_version = read_string("MEDAKA_VERSION")
    }
    runtime {
      docker: "~{docker}"
      cpu: cpu
      memory: "~{memory} GB"
      disks: "local-disk " + disk_size + " HDD"
      disk: disk_size + " GB"
      maxRetries: 1
      preemptible: 0
    }
}
