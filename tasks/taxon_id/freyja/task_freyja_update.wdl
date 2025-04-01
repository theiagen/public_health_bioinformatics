version 1.0

task freyja_update_refs {
  input {
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freyja:1.5.3"
    Int disk_size = 25
    Int memory = 10
    Int cpu = 1
  }
  meta {
    volatile: true
  }
  command <<<
  # Create updated reference files
  mkdir freyja_update_refs 
  freyja update --outdir freyja_update_refs
    
  echo "Freyja reference files created using the freyja update command; Freyja Docker Image utilized: ~{docker}. More information can be found at https://github.com/andersen-lab/Freyja" > freyja_update_refs/update_log.txt

  >>>
  output {
    File updated_barcodes = "freyja_update_refs/usher_barcodes.feather"
    File updated_lineages = "freyja_update_refs/curated_lineages.json"
    File update_log = "freyja_update_refs/update_log.txt"
  }
  runtime {
    memory: memory + " GB"
    cpu: cpu
    docker: "~{docker}"
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
  }
}
