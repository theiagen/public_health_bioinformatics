version 1.0

task meningotype {
  meta {
    description: "Serotyping of Neisseria meningitidis"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/meningotype:0.8.5--pyhdfd78af_0"
    Int? cpu = 2

    # Parameters
    # --finetype      perform porA and fetA fine typing (default=off)
    # --porB          perform porB sequence typing (NEIS2020) (default=off)
    # --bast          perform Bexsero antigen sequence typing (BAST) (default=off)
    # --mlst          perform MLST (default=off)
    # --all           perform MLST, porA, fetA, porB, BAST typing (default=off)
    Boolean finetype = false
    Boolean porB = false
    Boolean bast = false
    Boolean mlst = false
    Boolean all = false
  }
  command <<<
    echo $(meningotype --version 2>&1) | sed 's/^.*meningotype v//' | tee VERSION
    meningotype \
      ~{true="--finetype" false="" finetype} \
      ~{true="--porB" false="" porB} \
      ~{true="--bast" false="" bast} \
      ~{true="--mlst" false="" mlst} \
      ~{true="--all" false="" all} \
      --cpus ~{cpu} \
      ~{assembly} \
      > ~{samplename}.tsv
  >>>
  output {
    File meningotype_results = "~{samplename}.tsv"
    String meningotype_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
