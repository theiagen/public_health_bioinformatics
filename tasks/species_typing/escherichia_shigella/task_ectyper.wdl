version 1.0

task ectyper {
  meta {
    description: "In-silico prediction of Escherichia coli serotype"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/ectyper:1.0.0--pyhdfd78af_1"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 8

    # ECTyper Parameters
    #  --opid           [integer] Percent identity required for an O antigen allele match [default: 90]
    #  --opcov          [integer] Minumum percent coverage required for an O antigen allele match [default: 90]
    #  --hpid           [integer] Percent identity required for an H antigen allele match [default: 95]
    #  --hpcov          [integer] Minumum percent coverage required for an H antigen allele match [default: 50]
    #  --verify         [boolean] Enable E. coli species verification
    #  --print_alleles  [boolean] Prints the allele sequences if enabled as the final column
    Int o_min_percent_identity = 90
    Int h_min_percent_identity = 95
    Int o_min_percent_coverage = 90
    Int h_min_percent_coverage = 50
    Boolean verify = false
    Boolean print_alleles = false
  }
  command <<<
    echo $(ectyper --version 2>&1) | sed 's/.*ectyper //; s/ .*\$//' | tee VERSION
    ectyper \
      ~{'-opid ' + o_min_percent_identity} \
      ~{'-hpid ' + h_min_percent_identity} \
      ~{'-opcov ' + o_min_percent_coverage} \
      ~{'-hpcov ' + h_min_percent_coverage} \
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
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
