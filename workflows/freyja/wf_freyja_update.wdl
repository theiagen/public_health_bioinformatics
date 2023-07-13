version 1.0

workflow freyja_update {
  input {
    String gcp_uri
  }
  call freyja_update_refs {
    input:
  }
  call transfer_files {
    input:
      updated_barcodes = freyja_update_refs.updated_barcodes,
      updated_lineages = freyja_update_refs.updated_lineages,
      update_log = freyja_update_refs.update_log,
      gcp_uri = gcp_uri
  }
  output {
  }
}
task freyja_update_refs {
  input {
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:1.3.4"
    Int disk_size = 100
  }
  meta {
    volatile: true
  }
  command <<<
  # Create updated refrence files
  mkdir freyja_update_refs 
  freyja update --outdir freyja_update_refs
    
  echo "Freyja reference files created using the freyja update command; Freyja Docker Image utilized: ~{docker}. More information can be found at https://github.com/andersen-lab/Freyja" > freyja_update_refs/update_log.txt
   
  >>>
  runtime {
    memory: "16 GB"
    cpu: 4
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
  output {
    File updated_barcodes = "freyja_update_refs/usher_barcodes.csv"
    File updated_lineages = "freyja_update_refs/curated_lineages.json"
    File update_log = "freyja_update_refs/update_log.txt"
  }
}
task transfer_files {
  input {
    String gcp_uri
    File updated_barcodes
    File updated_lineages
    File update_log
    Int disk_size = 100
  }
  command <<<
  # transfer_files to specified gcp_uri
  date_tag=$(date +"%Y-%m-%d")
  
  gsutil -m cp ~{updated_barcodes} ~{updated_lineages} ~{update_log} ~{gcp_uri}/${date_tag}
  
  >>>
  runtime {
    memory: "4 GB"
    cpu: 2
    docker: "theiagen/utility:1.1"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
  output {    
  }
}