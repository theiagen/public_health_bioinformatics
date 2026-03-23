version 1.0

task augur_ancestral {
  input {
    File refined_tree
    File aligned_fasta
    String build_name

    String inference = "joint" # options: joint, marginal
    Boolean keep_ambiguous = false # do not infer nucleotides at ambiguous (N) sites
    Boolean infer_ambiguous = false # infer nucleotides at ambiguous sites and replace with most likely
    Boolean keep_overhangs = false # do not infer nucleotides for gaps on either side of the alignment

    Int disk_size = 50
    Int memory = 50
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/augur:31.5.0"
  }
  command <<<
    AUGUR_RECURSION_LIMIT=10000 augur ancestral \
      --tree "~{refined_tree}" \
      --alignment "~{aligned_fasta}" \
      --output-node-data "~{build_name}_nt_muts.json" \
      --output-sequences "~{build_name}_ancestral_sequences.fasta" \
      --inference ~{default="joint" inference} \
      ~{true="--keep-ambiguous" false="" keep_ambiguous} \
      ~{true="--infer-ambiguous" false="" infer_ambiguous} \
      ~{true="--keep-overhangs" false="" keep_overhangs} 
  >>>
  output {
    File ancestral_nt_muts_json = "~{build_name}_nt_muts.json"
    File ancestral_sequences = "~{build_name}_ancestral_sequences.fasta"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" 
    dx_instance_type: "mem3_ssd1_v2_x8"
    preemptible: 0
    maxRetries: 3
  }
}