version 1.0

task download_terra_table {
  meta {
    description: "This task downloads a Terra table and reduces it to only include the samples of interest."   
    
    # added so that call caching is always turned off
    volatile: true
  }
  input {
    String terra_table_name
    String terra_workspace_name
    String terra_project_name
    Int disk_size = 10
    Int memory = 1
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21"
  }
  command <<<
    python3 /scripts/export_large_tsv/export_large_tsv.py --project ~{terra_project_name} --workspace ~{terra_workspace_name} --entity_type ~{terra_table_name} --tsv_filename "~{terra_table_name}.tsv"
  >>>
  output {
    File terra_table = "~{terra_table_name}.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
  }
}