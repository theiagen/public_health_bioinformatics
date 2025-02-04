version 1.0

task krona {
  input {
    File kraken2_report
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/krona:2.8.1"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  command <<<
    # Get VERSION
    ktImportText | grep "KronaTools" | cut -d' ' -f3 | tee VERSION

    # run KrakenTools kreport2krona.py to enable viral compatibility
    kreport2krona.py -r ~{kraken2_report} -o ~{samplename}_krona.txt

    # Run krona with taxonomy on krakren report
    ktImportText ~{samplename}_krona.txt -o ~{samplename}_krona.html
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