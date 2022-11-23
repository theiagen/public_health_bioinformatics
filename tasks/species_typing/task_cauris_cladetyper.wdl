version 1.0

task cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
    Int kmer_size = 11
    String docker_image = "quay.io/biocontainers/hesslab-gambit:0.5.1--py37h8902056_0"
    Int memory = 16
    Int cpu = 8
    Int disk_size = 100
    File ref_clade1 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_reference.fasta"
    String ref_clade1_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_GCA_002759435_Cauris_B8441_V2_genomic.gbff"
    File ref_clade2 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_reference.fasta"
    String ref_clade2_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_CP043531.1.B11220.gb"
    File ref_clade3 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_reference.fasta"
    String ref_clade3_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.gbff"
    File ref_clade4 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_reference.fasta"
    String ref_clade4_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_PYGM01000135.1.B11243.gb"
    File ref_clade5 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_reference.fasta"
    String ref_clade5_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_CP050673.1.IFRC2087.gb"
  }
  command <<<
    # date and version control
    date | tee DATE
    gambit --version | tee VERSION

    gambit signatures create -o my-signatures.h5 -k ~{kmer_size} -p ATGAC ~{ref_clade1} ~{ref_clade2} ~{ref_clade3} ~{ref_clade4} ~{ref_clade5} ~{assembly_fasta}
    gambit dist --qs my-signatures.h5 --square -o ~{samplename}_matrix.csv
    
    cat ~{samplename}_matrix.csv | sort -k7 -t ',' | head -3 | tail -1 | rev | cut -d '/' -f1 | rev | cut -d ',' -f1 | cut -d '_' -f2 | tee CLADETYPE
    
    if [$CLADETYPE  
  >>>
  output {
    String gambit_cladetype = read_string("CLADETYPE")
    String clade_spec_ref = read_string("CLADEREF")
    String date = read_string("DATE")
    String version = read_string("VERSION")
    String gambit_cladetyper_docker_image = docker_image
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
