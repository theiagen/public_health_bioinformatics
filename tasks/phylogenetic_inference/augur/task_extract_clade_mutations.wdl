version 1.0

task extract_clade_mutations{
  input {
    File tree
    File metadata_tsv
    String clade_columns
    String tip_column
    File nt_mutations
    File? aa_mutations

    Boolean allow_missing_tips = true
    Boolean skip_singleton_clades = true

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/extract_clade_mutations:0.0.1"
    Int disk_size = 10
    Int memory = 4
    Int cpu = 1
  }
  String clades_file_path = basename(tree) + ".clades.tsv"
  command <<<
    # set -euo pipefail to avoid silent failure
    set -euo pipefail

    extract_nextclades.py \
      --tree ~{tree} \
      --metadata ~{metadata_tsv} \
      --clade_cols ~{clade_cols} \
      --tip_col ~{tip_col} \
      --output_tsv ~{clades_file_path}
      --nt_mutations ~{nt_mutations} \
      ~{if defined(aa_mutations) then "--aa_mutations " + aa_mutations else ""} \
      ~{true="--noncomprehensive" false="" allow_missing_tips} \
      ~{true="--skip_singletons" false="" skip_singleton_clades}
      -o ~{clades_file_path}
  >>>
  output {
    File clades_file = "~{clades_file_path}"
  }
  runtime {
    docker: docker 
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}
