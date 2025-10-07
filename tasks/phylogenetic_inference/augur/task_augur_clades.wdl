version 1.0

task augur_clades {
  input {
    File refined_tree
    File ancestral_nt_muts_json
    File translated_aa_muts_json
   # File reference_fasta
    String build_name

    File clades_tsv # tsv file containing clade definitions by amino acid
    Int disk_size = 50
    Int memory = 2
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur clades \
      --tree "~{refined_tree}" \
      --mutations "~{ancestral_nt_muts_json}" "~{translated_aa_muts_json}" \
      --clades "~{clades_tsv}" \
      --output-node-data "~{build_name}_clades.json"

    #--reference "{reference_fasta}" \
  >>>
  output {
    File clade_assignments_json = "~{build_name}_clades.json"
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