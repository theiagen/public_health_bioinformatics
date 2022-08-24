version 1.0

task ectyper {
  meta {
    description: "In-silico prediction of Escherichia coli serotype"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/ectyper:1.0.0--pyhdfd78af_1"
    Int? cpu = 4

    # ECTyper Parameters
    #  --opid           [integer] Percent identity required for an O antigen allele match [default: 90]
    #  --opcov          [integer] Minumum percent coverage required for an O antigen allele match [default: 90]
    #  --hpid           [integer] Percent identity required for an H antigen allele match [default: 95]
    #  --hpcov          [integer] Minumum percent coverage required for an H antigen allele match [default: 50]
    #  --verify         [boolean] Enable E. coli species verification
    #  --print_alleles  [boolean] Prints the allele sequences if enabled as the final column
    Int opid = 90
    Int hpid = 95
    Int opcov = 90
    Int hpcov = 50
    Boolean verify = false
    Boolean print_alleles = false
  }
  command <<<
    echo $(ectyper --version 2>&1) | sed 's/.*ectyper //; s/ .*\$//' | tee VERSION
    ectyper \
      ~{'-opid ' + opid} \
      ~{'-hpid ' + hpid} \
      ~{'-opcov ' + opcov} \
      ~{'-hpcov ' + hpcov} \
      ~{true="--verify" false="" verify} \
      ~{true="-s" false="" print_alleles} \
      --cores ~{cpu} \
      --output ./ \
      --input ~{assembly}
    mv output.tsv ~{samplename}.tsv
    # parse ECTyper TSV
    cut -f 5 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_SEROTYPE
  >>>
  output {
    File ectyper_results = "~{samplename}.tsv"
    String ectyper_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String ectyper_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
