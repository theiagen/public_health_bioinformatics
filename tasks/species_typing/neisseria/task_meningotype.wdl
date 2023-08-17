version 1.0

task meningotype {
  meta {
    description: "Serotyping of Neisseria meningitidis"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/meningotype:0.8.5--pyhdfd78af_0"
    Int disk_size = 100
    Int cpu = 2
  }
  command <<<
    # get version information
    echo $(meningotype --version 2>&1) | sed 's/^.*meningotype v//' | tee VERSION

    # Parameters
    # --finetype      perform porA and fetA fine typing (default=off)
    # --porB          perform porB sequence typing (NEIS2020) (default=off)
    # --bast          perform Bexsero antigen sequence typing (BAST) (default=off)
    # --mlst          perform MLST (default=off)
    # --all           perform MLST, porA, fetA, porB, BAST typing (default=off)

    meningotype \
      --finetype \
      --porB \
      --bast \
      --cpus ~{cpu} \
      ~{assembly} \
      > ~{samplename}.tsv

    tail -1 ~{samplename}.tsv | awk '{print $2}' | tee MENINGOTYPE_SEROTYPE
    tail -1 ~{samplename}.tsv | awk '{print $5}' | tee MENINGOTYPE_PORA
    tail -1 ~{samplename}.tsv | awk '{print $6}' | tee MENINGOTYPE_FETA
    tail -1 ~{samplename}.tsv | awk '{print $7}' | tee MENINGOTYPE_PORB
    tail -1 ~{samplename}.tsv | awk '{print $8}' | tee MENINGOTYPE_FHBP
    tail -1 ~{samplename}.tsv | awk '{print $9}' | tee MENINGOTYPE_NHBA
    tail -1 ~{samplename}.tsv | awk '{print $10}' | tee MENINGOTYPE_NADA
    tail -1 ~{samplename}.tsv | awk '{print $11}' | tee MENINGOTYPE_BAST

  >>>
  output {
    File meningotype_tsv = "~{samplename}.tsv"
    String meningotype_version = read_string("VERSION")
    String meningotype_serogroup = read_string("MENINGOTYPE_SEROTYPE")
    String meningotype_PorA = read_string("MENINGOTYPE_PORA")
    String meningotype_FetA = read_string("MENINGOTYPE_FETA")
    String meningotype_PorB = read_string("MENINGOTYPE_PORB")
    String meningotype_fHbp = read_string("MENINGOTYPE_FHBP")
    String meningotype_NHBA = read_string("MENINGOTYPE_NHBA")
    String meningotype_NadA = read_string("MENINGOTYPE_NADA")
    String meningotype_BAST = read_string("MENINGOTYPE_BAST")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
  }
}
