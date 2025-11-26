version 1.0

task qualimap {
  input {
    String samplename
    File bam_file
    Int disk_size = 50
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/qualimap-custom-html:2.3"
  }

  command <<<
    set -euo pipefail

    # set tmp dir for jaba
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp

    # get the version
    tee VERSION
    qualimap --help | grep "QualiMap v" >> VERSION
    # run qualimap bamqc for viz
    qualimap bamqc -bam ~{bam_file} -outdir qualimap_results

    # extract raw data for custom scripts/plots
    cp qualimap_results/raw_data_qualimapReport/coverage_across_reference.txt ~{samplename}_genome_coverage_across_reference.txt
    cp qualimap_results/raw_data_qualimapReport/mapping_quality_across_reference.txt ~{samplename}_mapping_quality_across_reference.txt
    # create my custom interactive plots
    plot_coverage.py ~{samplename}_genome_coverage_across_reference.txt ~{samplename}_mapping_quality_across_reference.txt ~{samplename}
    # zip results
    tar -zcvf ~{samplename}_qualimap_reports.tar.gz qualimap_results
  >>>
  output {
    String version = read_string("VERSION")
    String qualimap_docker = docker
    File qualimap_reports_bundle = "~{samplename}_qualimap_reports.tar.gz"
    File qualimap_coverage_plots_html = "~{samplename}_coverage_plot.html"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
    maxRetries: 3
  }
}