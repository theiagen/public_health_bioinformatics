version 1.0

task racon {
  input {
    File unpolished_fasta
    File read1
    String samplename
    Int polishing_rounds = 1 # Default: 1 polishing round
    Int cpu = 8                  
    Int memory = 32            
    Int disk_size = 100 
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/racon:1.5.0-minimap2-generic-v3"
  }
  command <<< 
    set -euo pipefail

    minimap2 --version | tee MINIMAP2_VERSION
    racon --version | tee -a RACON_VERSION

    echo "Starting Racon polishing process..."

    # Initialize the input assembly
    intermediate_fasta="~{unpolished_fasta}"

    # Loop through polishing rounds
    for i in $(seq 1 ~{polishing_rounds}); do
      echo "Polishing round $i..."

      # Align reads to the current assembly
      minimap2 -x map-ont -t ~{cpu} "${intermediate_fasta}" "~{read1}" > "~{samplename}_round${i}.paf"

      # Run Racon to polish the assembly
      racon \
        -t ~{cpu} \
        "~{read1}" \
        "~{samplename}_round${i}.paf" \
        "${intermediate_fasta}" \
        > "~{samplename}_round${i}.polished.fasta"

      # Update for the next round
      intermediate_fasta="~{samplename}_round${i}.polished.fasta"
    done

    # Save the final polished assembly
    mv "${intermediate_fasta}" "~{samplename}_final_polished.fasta"
    echo "Polishing complete. Final assembly saved as ~{samplename}_final_polished.fasta"
  >>>
  output {
    File polished_fasta = "~{samplename}_final_polished.fasta"  
    String racon_version = read_string("RACON_VERSION")   
    String minimap2_version = read_string("MINIMAP2_VERSION") 
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
    preemptible: 0
  }
}

