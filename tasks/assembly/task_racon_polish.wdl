version 1.0

task racon_polish {
    input {
        File assembly_fasta 
        File read1 
        Int racon_rounds = 1
        Int cpu = 4                  
        Int memory = 16                   
        Int disk_size = 100                     
        String samplename                  
        String minimap2_docker = "staphb/minimap2:2.24" 
        String racon_docker = "staphb/racon:1.4.20"      
    }

    command <<< 
        # Initialize variables
        cp ~{assembly_fasta} ~{samplename}.unpolished.fasta
        POLISHED_FASTA=~{samplename}.unpolished.fasta

        # Run Racon for the specified number of rounds
        for i in $(seq 1 ~{racon_rounds}); do
            echo "Polishing with Racon - Round $i of ~{racon_rounds}"
            mkdir -p racon_round_${i}

            # Step 1: Generate PAF alignments with Minimap2
            minimap2 -t ~{cpu} -x map-ont $POLISHED_FASTA ~{read1} \
                > racon_round_${i}/alignments.paf

            # Step 2: Run Racon using the generated PAF
            racon -t ~{cpu} ~{read1} racon_round_${i}/alignments.paf $POLISHED_FASTA \
                > racon_round_${i}/consensus.fasta

            # Update the polished FASTA for the next round
            POLISHED_FASTA=racon_round_${i}/consensus.fasta
        done

        # Copy final polished FASTA to the home directory with sample name
        cp $POLISHED_FASTA ~{samplename}.polished.fasta
    >>>

    output {
        File raconpolished_fasta = "~{samplename}.polished.fasta" 
    }

    runtime {
        docker: "~{racon_docker}" 
        cpu: cpu 
        memory: "~{memory} GB"  
        disks: "local-disk " + disk_size + " HDD" 
        disk: disk_size + " GB"                
        maxRetries: 3                         
        preemptible: 0                          
    }
}