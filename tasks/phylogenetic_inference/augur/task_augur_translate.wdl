version 1.0

task augur_translate {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    File reference_genbank
    String build_name

    File? genes # a file containing list of genes to translate (from nucleotides to amino acids)

    Int disk_size = 50
    Int mem_size = 32
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur translate \
      --tree "~{refined_tree}" \
      --ancestral-sequences "~{ancestral_nt_muts_json}" \
      --reference-sequence "~{reference_genbank}" \
      ~{"--genes " + genes} \
      --output-node-data "~{build_name}_aa_muts.json"
  >>>
  output {
    File translated_aa_muts_json = "~{build_name}_aa_muts.json"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/biocontainers/augur:22.0.2--pyhdfd78af_0"
    memory: mem_size + " GB"
    cpu : 1
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 0
    maxRetries: 3
  }
}