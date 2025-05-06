version 1.0

task quast {
  input {
    File assembly
    String samplename
    File? reference
    File? reference_gff
    Int min_contig_length = 500
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/quast:5.0.2"
    Int disk_size = 100
    Int memory = 2 # added default value
    Int cpu = 2 # added default value
  }
  command <<<
    # capture date and version
    date | tee DATE
    quast.py --version | grep QUAST | tee VERSION

    quast.py \
      ~{assembly} \
      -o . \
      --min-contig ~{min_contig_length} \
      ~{'-r ' + reference} \
      ~{'-g ' + reference_gff}

    mv report.tsv ~{samplename}_report.tsv

    # Check if required files/directories for Icarus exist
    if [[ -d "icarus_viewers" ]]; then
      # Create tgz of files needed for Icarus HTML viewer
      mkdir ~{samplename}_icarus_viewer
      mv icarus_viewers ~{samplename}_icarus_viewer/ && \
        mv icarus.html ~{samplename}_icarus_viewer/ && 
        mv report.html ~{samplename}_icarus_viewer/
      tar -zcvf ~{samplename}_icarus_viewer.tgz ~{samplename}_icarus_viewer
    fi

    python <<CODE
    import csv
    #grab output genome length and number contigs by column header
    with open("~{samplename}_report.tsv",'r') as tsv_file:
      tsv_reader = csv.reader(tsv_file, delimiter="\t")
      for line in tsv_reader:
          if "Total length" in line[0]:
            with open("GENOME_LENGTH", 'wt') as genome_length:
              genome_length.write(line[1])
          if "# contigs" in line[0]:
            with open("NUMBER_CONTIGS", 'wt') as number_contigs:
              number_contigs.write(line[1])
          if "N50" in line[0]:
            with open("N50_VALUE", 'wt') as n50_value:
              n50_value.write(line[1])
          if "GC" in line[0]:
            with open("GC_PERCENT", 'wt') as gc_percent:
              gc_percent.write(line[1])
          if "Largest contig" in line[0]:
            with open("LARGEST_CONTIG", 'wt') as largest_contig:
              largest_contig.write(line[1])
          if "# N's per 100 kbp" in line[0]:
            with open("UNCALLED_BASES", "wt") as uncalled_bases:
              uncalled_bases.write(line[1])
    CODE

  >>>
  output {
    File quast_report = "~{samplename}_report.tsv"
    File? icarus_report = "~{samplename}_icarus_viewer.tgz"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    Int genome_length = read_int("GENOME_LENGTH")
    Int number_contigs = read_int("NUMBER_CONTIGS")
    Int n50_value = read_int("N50_VALUE")
    Float gc_percent = read_float("GC_PERCENT")
    Int largest_contig = read_int("LARGEST_CONTIG")
    Float uncalled_bases = read_float("UNCALLED_BASES")
    String quast_docker = docker    
  }
  runtime {
    docker:  "~{docker}"
    memory:  "~{memory} GB"
    cpu:   "~{cpu}"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}