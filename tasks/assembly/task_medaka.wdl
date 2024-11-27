version 1.0

task medaka_consensus {
    input {
      File assembly_fasta
      String samplename
      File read1
      String medaka_model
      Int polish_rounds = 1
      Int cpu = 4
      Int memory = 16
      Int disk_size = 100
      String docker = "us-docker.pkg.dev/general-theiagen/staphb/medaka:2.0.1"
    }
    command <<< 
      medaka --version | tee VERSION

      # Initialize the input for polishing with the provided assembly FASTA
        cp ~{assembly_fasta} polished_input.fasta
        cp ~{assembly_fasta} polished_input.fasta

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
                --threads ~{cpu}

            # Update the input for the next round
            cp "$polish_dir/consensus.fasta" polished_input.fasta
        done

        # Rename final polished outputs with sample name
        mv polished_input.fasta ~{samplename}.polished.fasta
        mv polish_round_~{polish_rounds}/consensus.vcf.gz ~{samplename}.polished.vcf.gz
    >>>
    output {
      File medaka_fasta = "~{samplename}.polished.fasta"
      File medaka_vcf = "~{samplename}.polished.vcf.gz"
      String medaka_version = read_string("VERSION")
    }
    runtime {
      docker: "~{docker}"
      cpu: cpu
      memory: "~{memory} GB"
      disks: "local-disk " + disk_size + " HDD"
      disk: disk_size + " GB"
      maxRetries: 3
      preemptible: 0
    }
}
