version 1.0

task sistr {
  meta {
    description: "Serovar prediction of Salmonella assemblies"
  }
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/sistr:1.1.3"
    Int cpu = 2
    Int memory = 8
    Int disk_size = 100

    # Parameters
    # --use-full-cgmlst-db  Use the full set of cgMLST alleles which can include highly similar alleles. By default the smaller "centroid" alleles or representative alleles are used for each marker. 
    Boolean use_full_cgmlst_db = false
  }
  command <<<
    set -euo pipefail

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
    
    # rename to .tsv suffix
    mv ~{samplename}.tab ~{samplename}.tsv
    
    # parse sistr TSV
    cut -f 16 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_SEROTYPE
    cut -f 15 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_SEROGROUP
    cut -f 10 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_H1_ANTIGEN
    cut -f 11 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_H2_ANTIGEN
    cut -f 12 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_O_ANTIGEN
    cut -f 18 ~{samplename}.tsv | tail -n 1 | tee PREDICTED_SEROTYPE_CGMLST
    cut -f 1 ~{samplename}.tsv | tail -n 1 | tee ANTIGENIC_FORMULA
    
  >>>
  output {
    File sistr_results = "~{samplename}.tsv"
    File sistr_allele_json = "~{samplename}-allele.json"
    File sistr_allele_fasta = "~{samplename}-allele.fasta"
    File sistr_cgmlst = "~{samplename}-cgmlst.csv"
    String sistr_version = read_string("VERSION")
    String sistr_antigenic_formula = read_string("ANTIGENIC_FORMULA")
    String sistr_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String sistr_serogroup = read_string("PREDICTED_SEROGROUP")
    String sistr_h1_antigens = read_string("PREDICTED_H1_ANTIGEN")
    String sistr_h2_antigens = read_string("PREDICTED_H2_ANTIGEN")
    String sistr_o_antigens = read_string("PREDICTED_O_ANTIGEN")
    String sistr_serotype_cgmlst = read_string("PREDICTED_SEROTYPE_CGMLST")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
