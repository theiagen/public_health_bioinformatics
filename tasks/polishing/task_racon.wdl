version 1.0

task racon {
  input {
    File unpolished_fasta
    File read1
    Int polishing_rounds = 1      # Default: 1 polishing round
    Int cpu = 8                  
    Int memory = 16             
    Int disk_size = 100   
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/racon:1.5.0-minimap2"
  }
  command <<< 
    set -euo pipefail

    minimap2 --version | tee MINIMAP2_VERSION
    racon --version | tee -a RACON_VERSION

    # Start with the unpolished FASTA as the input
    intermediate_fasta="~{unpolished_fasta}"

     for i in $(seq 1 ~{polishing_rounds}); do
       echo "Starting Medaka polishing round $i"

      # Generate alignments with Minimap2
      minimap2 -t ~{cpu} "${intermediate_fasta}" "~{read1}" > "~{samplename}_round${i}.paf"

      # Run Racon for polishing
      racon \
        -t ~{cpu} \
        "~{read1}" \
        "~{samplename}_round${i}.paf" \
        "${intermediate_fasta}" \
        > "~{samplename}_round${i}.polished.fasta"

      # Update current_fasta for the next round
      intermediate_fasta="~{samplename}_round${i}.polished.fasta"
    done

    # Move the final polished assembly to the output
    mv "${intermediate_fasta}" "~{samplename}_final_polished.fasta"
  >>>
  output {
    File polished_fasta = "~{samplename}_final_polished.fasta"  
    String racon_version = read_string("RACON_VERSION")   
    String minimap2_version = read_string("MINIMAP2_VERSION") 
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 3
    preemptible: 0
  }
}

