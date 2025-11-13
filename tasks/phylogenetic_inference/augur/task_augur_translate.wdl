version 1.0

task augur_translate {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    File reference_genbank
    String build_name

    File? genes # a file containing list of genes to translate (from nucleotides to amino acids)

    Int disk_size = 50
    Int cpu = 1
    Int memory = 32
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
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
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    preemptible: 0
    maxRetries: 3
  }
}