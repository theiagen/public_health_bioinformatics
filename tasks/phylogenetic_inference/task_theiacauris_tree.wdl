version 1.0

task theiacauris_mashtree_fasta {
  input {
    File assembly_fasta
    String cluster_name
    Int truncLength = 250
    String sort_order = "ABC"
    Int genomesize = 5000000
    Int mindepth = 5
    Int kmerlength = 21
    Int sketchsize = 10000
    Int cpu = 8
    Int memory = 32
    File ref_clade1 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade_I_reference.fasta"
    File ref_clade2 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade_II_reference.fasta"
    File ref_clade3 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade_III_reference.fasta"
    File ref_clade4 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade_IV_reference.fasta"
    File ref_clade5 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade_V_reference.fasta"
    File ref_other = "gs://theiagen-public-files/terra/candida_auris_refs/candida_auris_B11221_PGLS01000001.1.fasta"
  }
  command <<<
    # date and version control
    date | tee DATE
    mashtree -v | tee VERSION
    
    # organize input assemblies
    mkdir mash_assemblies
    if [[ ! -z ~{ref_other} ]]; then
      mv ~{ref_other} mash_assemblies
    fi
    mv ~{ref_clade1} mash_assemblies
    mv ~{ref_clade2} mash_assemblies
    mv ~{ref_clade3} mash_assemblies
    mv ~{ref_clade4} mash_assemblies
    mv ~{ref_clade5} mash_assemblies
    mv ~{assembly_fasta} mash_assemblies
    #run mashtree
    mashtree \
      ~{'--truncLength ' + truncLength} \
      ~{'--sort-order ' + sort_order} \
      ~{'--genomesize ' + genomesize} \
      ~{'--mindepth ' + mindepth} \
      ~{'--kmerlength ' + kmerlength} \
      ~{'--sketch-size ' + sketchsize} \
      ~{'--numcpus ' + cpu} \
      ~{'--outmatrix ' + cluster_name + '.tsv'} \
      ~{'--outtree ' + cluster_name + '.nwk'} \
      mash_assemblies/*
      
      #Find and return min value col header of the min dist 
      
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File mashtree_matrix = "~{cluster_name}.tsv"
    File mashtree_tree = "~{cluster_name}.nwk"
  }
  runtime {
    docker: "quay.io/staphb/mashtree:1.2.0"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
