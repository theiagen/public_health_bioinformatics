version 1.0

task cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
    Int kmer_size = 11
    String docker_image = "quay.io/biocontainers/hesslab-gambit:0.5.1--py37h8902056_0"
    Int memory = 16
    Int cpu = 8
    File ref_clade1 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_GCA_002759435.2_Cand_auris_B8441_V2_genomic.fna"
    String ref_clade1_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade1_GCA_002759435_Cauris_B8441_V2_genomic.gbff"
    File ref_clade2 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_reference.fasta"
    String ref_clade2_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade2_CP043531.1.B11220.gb"
    File ref_clade3 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_reference.fasta"
    String ref_clade3_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.gbff"
    File ref_clade4 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_reference.fasta"
    String ref_clade4_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.gbff"
    File ref_clade5 = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.fasta"
    String ref_clade5_annotated = "gs://theiagen-public-files/terra/candida_auris_refs/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.gbff"
  }
  command <<<
    # date and version control
    date | tee DATE
    gambit --version | tee VERSION
    # create gambit signature file for five clades + input assembly
    gambit signatures create -o my-signatures.h5 -k ~{kmer_size} -p ATGAC ~{ref_clade1} ~{ref_clade2} ~{ref_clade3} ~{ref_clade4} ~{ref_clade5} ~{assembly_fasta}
    # calculate distance matrix for all six signatures
    gambit dist --qs my-signatures.h5 --square -o ~{samplename}_matrix.csv
    # parse matrix to see closest clade to input assembly
    ## sort by 7th column (distance against input sequence)
    ## take top three columns: header, header, top hit
    ## take bottom of these three rows (top hit)
    ## grab only file name
    top_clade=$(sort -k7 -t ',' "~{samplename}_matrix.csv" | head -3 | tail -n-1 | awk -F',' '{print$1}')
    #clade_type=$(cat CLADETYPE)
    if [ "${top_clade}" == "~{ref_clade1}" ] ; then
      echo "~{ref_clade1_annotated}" > CLADEREF
      echo "Clade1" > CLADETYPE
    elif [ "${top_clade}" == "~{ref_clade2}" ] ; then
      echo "~{ref_clade2_annotated}" > CLADEREF
      echo "Clade2" > CLADETYPE
    elif [ "${top_clade}" == "~{ref_clade3}" ] ; then
      echo "~{ref_clade3_annotated}" > CLADEREF
      echo "Clade3" > CLADETYPE
    elif [ "${top_clade}" == "~{ref_clade4}" ] ; then
      echo "~{ref_clade4_annotated}" > CLADEREF
      echo "Clade4" > CLADETYPE
    elif [ "${top_clade}" == "~{ref_clade5}" ] ; then
      echo "~{ref_clade5_annotated}" > CLADEREF
      echo "Clade5" > CLADETYPE
    else
      echo "None" > CLADEREF
    fi

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
