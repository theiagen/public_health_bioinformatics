version 1.0

task staphopiasccmec {
  meta {
    description: "Primer based SCCmec typing of Staphylococcus aureus genomes"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/staphopia-sccmec:1.0.0--hdfd78af_0"
    Int? cpu = 2

    # Parameters
    # --hamming Report the results as hamming distances
    Boolean hamming = false
  }
  command <<<
    staphopia-sccmec --version 2>&1 | sed 's/^.*staphopia-sccmec //' | tee VERSION
    staphopia-sccmec \
      ~{true="--hamming" false="" hamming} \
      --assembly ~{assembly} > ~{samplename}.tsv
  >>>
  output {
    File staphopiasccmec_results = "~{samplename}.tsv"
    String staphopiasccmec_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
