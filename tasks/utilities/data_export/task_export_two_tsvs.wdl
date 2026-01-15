version 1.0

task export_two_tsvs {
  input {
    String terra_project1
    String? terra_project2
    String terra_workspace1
    String? terra_workspace2
    String datatable1
    String datatable2
    Int disk_size = 10
    Int memory = 1
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    set -euo pipefail
    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project1} --workspace ~{terra_workspace1} --entity_type ~{datatable1} --tsv_filename "~{datatable1}_table1.tsv"

    # check if second project is provided; if not, use first
    if [[ -z "~{terra_project2}" ]]; then
      PROJECT2="~{terra_project1}"
    else
      PROJECT2="~{terra_project2}"
    fi

    # check if second workspace is provided; if not, use first
    if [[ -z "~{terra_workspace2}" ]]; then
      WORKSPACE2="~{terra_workspace1}"
    else
      WORKSPACE2="~{terra_workspace2}"
    fi

    python3 /scripts/export_large_tsv/export_large_tsv.py --project ${PROJECT2} --workspace ${WORKSPACE2} --entity_type ~{datatable2} --tsv_filename "~{datatable2}_table2.tsv"

    if [[ $(wc -l ~{datatable1}_table1.tsv | cut -f1 -d' ') -eq $(wc -l ~{datatable2}_table2.tsv | cut -f1 -d' ') ]]; then
      echo true | tee CONTINUE
    else 
      echo false | tee CONTINUE
    fi
  >>>
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
  output {
    File datatable1_tsv = "~{datatable1}_table1.tsv"
    File datatable2_tsv = "~{datatable2}_table2.tsv"
    Boolean same_table_length = read_boolean("CONTINUE")
  }
}
