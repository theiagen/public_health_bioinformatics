version 1.0

task krona {
  input {
    File kraken2_report
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/krona:2.7.1--pl526_5"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
    # Get VERSION
    ktImportTaxonomy 2>&1 | sed -n '/KronaTools /p' | sed 's/^.*KronaTools //; s/ - ktImportTaxonomy.*//' | tee VERSION

    # Get taxonomy file 
    ktUpdateTaxonomy.sh taxonomy

    # Run krona with taxonomy on kraken2 report
    ktImportTaxonomy -o ~{samplename}_krona.html ~{kraken2_report} -tax taxonomy -t 5
  >>>
  output {
    String krona_version = read_string("VERSION")
    String krona_docker = docker
    File krona_html = "~{samplename}_krona.html"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      disk: disk_size + " GB"
      preemptible: 0
  }
}