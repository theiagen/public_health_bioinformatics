version 1.0

# Notes:
# Gamma is run 3 times within OHs workflow, one for ar, one for hv, and one for plasmid replicons.
# Run gamma 3 times, HV, AR, and plasmid replicons by passing the appropriate database
# and call gamma as gamma_hv, gamma_ar, and gamma_rep in the workflow.
# 
# GAMMA_S runs on nucleotides and does not translate to protein, like normal gamma, 
# this can be chosen with run_gammas. default of false.

task gamma {
  input {
    File assembly
    File gamma_hv_db
    File gamma_ar_db
    File gamma_rep_db
    Int min_percent_identity = 90
    Boolean run_gammas = false
    Boolean output_gff = true
    Boolean output_fasta = true
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/gamma:2.2"
    Int cpu = 2
    Int disk_size = 50
    Int memory = 8
  }
  command <<<
    if ~{run_gammas}; then
      # Hyper Virulence Genes
      if ~{gamma_hv_db}; then
        GAMMA_S.py ~{assembly} \
          ~{gamma_hv_db} \
          ~{samplename + "_hv"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name 
      fi
      # AR genes
      if ~{gamma_ar_db}; then
        GAMMA_S.py ~{assembly} \
          ~{gamma_ar_db} \
          ~{samplename + "_ar"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name
      fi
      # Plasmid replicons
      if ~{gamma_rep_db}; then
        GAMMA_S.py ~{assembly} \
          ~{gamma_rep_db} \
          ~{samplename + "_rep"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name
      fi
    else
      # Hyper Virulence Genes
      if ~{gamma_hv_db}; then
        GAMMA.py ~{assembly} \
          ~{gamma_hv_db} \
          ~{samplename + "_hv"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name 
      fi
      # AR Genes
      if ~{gamma_ar_db}; then
        GAMMA.py ~{assembly} \
          ~{gamma_ar_db} \
          ~{samplename + "_ar"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name
      fi
      # Plasmid replicons
      if ~{gamma_rep_db}; then
        GAMMA.py ~{assembly} \
          ~{gamma_rep_db} \
          ~{samplename + "_rep"} \
          $(if ~{output_gff}; then echo "--gff"; fi) \
          $(if ~{output_fasta}; then echo "--fasta"; fi) \
          --identity ~{min_percent_identity} \
          --name
      fi
    fi
  >>>
  output {
    String? gamma_hv_results = "~{samplename}_hv.gamma"
    String? gamma_ar_results = "~{samplename}_ar.gamma"
    String? gamma_rep_results = "~{samplename}_rep.gamma"
    String? gamma_hv_gff = "~{samplename}_hv.gff"
    String? gamma_ar_gff = "~{samplename}_ar.gff"
    String? gamma_rep_gff = "~{samplename}_rep.gff"
    String? gamma_hv_fasta = "~{samplename}_hv.fasta"
    String? gamma_ar_fasta = "~{samplename}_ar.fasta"
    String? gamma_rep_fasta = "~{samplename}_rep.fasta"
    String gamma_docker = docker
  }
  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}