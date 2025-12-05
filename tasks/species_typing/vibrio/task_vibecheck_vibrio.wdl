version 1.0

task vibecheck_vibrio {
  meta {
    description: "Rapidly assigns O1 Vibrio cholerae sequencing data to transmission lineages using variant frequency demixing."
  }
  input {
    File read1
    File? read2
    File? lineage_barcodes
    Float subsampling_fraction = 0.2
    Boolean skip_subsampling = false
    String docker = "watronfire/vibecheck:2025.02.24"
    Int disk_size = 16
    Int cpu = 2
    Int memory = 3
  }
  command <<<
    set -ex
   
    # Capture Vibecheck version information.
    vibecheck -v | tee VERSION

    # Unclear if this will work with single-reads, will explore and update.
    vibecheck ~{read1} ~{read2} \
        --outdir . \
        ~{"--barcodes " + lineage_barcodes} \
        ~{"--subsampling_fraction " + subsampling_fraction} \
        ~{true='--no-detect' false='' skip_subsampling}

    # Parse Vibecheck results to capture lineage estimate, confidence, and any classification notes produced by Freyja.
    python3 <<CODE
    import csv
    with open( "lineage_report.csv", "rt" ) as csv_file:
      reader = csv.DictReader( csv_file )
      line = next( reader )
      for key in ["lineage", "confidence", "classification_notes"]:
        with open( key.upper(), "wt" ) as outf:
          outf.write(line[key])
    CODE
  >>>
  output {
    File vibecheck_lineage_report = "lineage_report.csv"
    String vibecheck_top_lineage = read_string("LINEAGE")
    Float vibecheck_confidence = read_float("CONFIDENCE")
    String vibecheck_classification_notes = read_string("CLASSIFICATION_NOTES")
    String vibecheck_version = read_string("VERSION")
    String vibecheck_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
