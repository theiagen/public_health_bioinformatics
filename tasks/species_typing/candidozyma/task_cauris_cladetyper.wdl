version 1.0

task cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
    Int kmer_size = 11
    
    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0"
    Int memory = 16

    File ref_clade1 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1_GCA_002759435.2_Cand_auris_B8441_V2_genomic.fasta"
    String ref_clade1_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1_GCA_002759435_Cauris_B8441_V2_genomic.gbff"
    File ref_clade2 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.fasta"
    String ref_clade2_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.gbff"
    File ref_clade3 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.fasta"
    String ref_clade3_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.gbff"
    File ref_clade4 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.fasta"
    String ref_clade4_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.gbff"
    File ref_clade5 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.fasta"
    String ref_clade5_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.gbff"
    File ref_clade6 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade6_GCA_032714025.1_ASM3271402v1_genomic.fasta"
    String? ref_clade6_annotated
    }
  command <<<
    gambit --version | tee VERSION

    # create gambit signature file for six clades + input assembly
    gambit signatures create -o my-signatures.h5 -k ~{kmer_size} -p ATGAC ~{ref_clade1} ~{ref_clade2} ~{ref_clade3} ~{ref_clade4} ~{ref_clade5} ~{ref_clade6} ~{assembly_fasta}
    # calculate distance matrix for all seven signatures
    gambit dist --qs my-signatures.h5 --square -o ~{samplename}_matrix.csv

    # parse matrix to see closest clade to input assembly
    ## sort by 8th column (distance against input sequence)
    ## take top three columns: header, header, top hit
    ## take bottom of these three rows (top hit)
    ## grab only file name (top_clade)
    top_clade=$(sort -k8 -t ',' "~{samplename}_matrix.csv" | head -3 | tail -n-1 | awk -F',' '{print$1}')

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
    elif [ "${top_clade}" == "~{ref_clade6}" ] ; then
      # clade 6 may not be defined
      if [ -z "~{ref_clade6_annotated}" ] ; then
        echo "None" > CLADEREF
      else
        echo "~{ref_clade6_annotated}" > CLADEREF
      fi
      echo "Clade6" > CLADETYPE
    else
      echo "None" > CLADEREF
      echo "" > CLADETYPE
    fi

  >>>
  output {
    String gambit_version = read_string("VERSION")
    String gambit_cladetype = read_string("CLADETYPE")
    String annotated_reference = read_string("CLADEREF")
    String gambit_cladetyper_docker_image = docker
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu    
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
