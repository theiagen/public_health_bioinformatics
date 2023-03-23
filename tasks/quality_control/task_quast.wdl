version 1.0

task quast {
  input {
    File assembly
    String samplename
    String docker = "quay.io/staphb/quast:5.0.2"
    Int disk_size = 100
    Int memory = 2 # added default value
    Int cpu = 2 # added default value
  }
  command <<<
    # capture date and version
    date | tee DATE
    quast.py --version | grep QUAST | tee VERSION

    quast.py ~{assembly} -o .
    mv report.tsv ~{samplename}_report.tsv
    
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

    CODE

  >>>
  output {
    File quast_report = "${samplename}_report.tsv"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    Int genome_length = read_int("GENOME_LENGTH")
    Int number_contigs = read_int("NUMBER_CONTIGS")
    Int n50_value = read_int("N50_VALUE")
    Float gc_percent = read_float("GC_PERCENT")    
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