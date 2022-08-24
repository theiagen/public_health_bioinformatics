version 1.0

task sistr {
  meta {
    description: "Serovar prediction of Salmonella assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/sistr_cmd:1.1.1--pyh864c0ab_2"
    Int? cpu = 4

    # Parameters
    # --use-full-cgmlst-db  Use the full set of cgMLST alleles which can include highly similar alleles. By default the smaller "centroid" alleles or representative alleles are used for each marker. 
    Boolean use_full_cgmlst_db = false
  }
  command <<<
    echo $(sistr --version 2>&1) | sed 's/^.*sistr_cmd //; s/ .*\$//' | tee VERSION
    sistr \
      --qc \
      ~{true="--use-full-cgmlst-db" false="" use_full_cgmlst_db} \
      --threads ~{cpu} \
      --alleles-output ~{samplename}-allele.json \
      --novel-alleles ~{samplename}-allele.fasta \
      --cgmlst-profiles ~{samplename}-cgmlst.csv \
      --output-prediction ~{samplename} \
      --output-format tab \
      ~{assembly}
    
    mv ~{samplename}.tab ~{samplename}.tsv
    
    # parse sistr TSV
    cut -f 15 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_SEROTYPE
    
  >>>
  output {
    File sistr_results = "~{samplename}.tsv"
    File sistr_allele_json = "~{samplename}-allele.json"
    File sistr_allele_fasta = "~{samplename}-allele.fasta"
    File sistr_cgmlst = "~{samplename}-cgmlst.csv"
    String sistr_version = read_string("VERSION")
    String sistr_predicted_serotype = read_string("PREDICTED_SEROTYPE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
