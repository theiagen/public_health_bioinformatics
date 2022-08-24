version 1.0

task spatyper {
  meta {
    description: "Computational method for finding spa types in Staphylococcus aureus"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/spatyper:0.3.3--pyhdfd78af_3"
    Int? cpu = 4

    # Parameters
    # --do_enrich Do PCR product enrichment
    Boolean do_enrich = false
  }
  command <<<
    echo \$(spaTyper --version 2>&1) | sed 's/^.*spaTyper //' | tee VERSION
    spaTyper \
      ~{true="--do_enrich" false="" do_enrich} \
      --fasta ~{assembly} \
      --output ~{samplename}.tsv
  >>>
  output {
      File spatyper_results = "~{samplename}.tsv"
      String spatyper_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
