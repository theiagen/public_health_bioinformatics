version 1.0

task fastqc_pe {
  input {
    File read1
    File read2
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String read2_name = basename(basename(basename(read2, ".gz"), ".fastq"), ".fq")
    Int? cpus = 2
    String docker="quay.io/staphb/fastqc:0.11.9"
  }
  command <<<
    # capture date and version
    date | tee DATE
    fastqc --version | grep FastQC | tee VERSION

    fastqc --outdir $PWD --threads ~{cpus} ~{read1} ~{read2}

    unzip -p ~{read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
    unzip -p ~{read2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

    READ1_SEQS=$(unzip -p ~{read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
    READ2_SEQS=$(unzip -p ~{read2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

    if [ $READ1_SEQS == $READ2_SEQS ]; then
      read_pairs=$READ1_SEQS
    else
      read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
    fi
    echo $read_pairs | tee READ_PAIRS
  >>>
  output {
    File fastqc1_html = "~{read1_name}_fastqc.html"
    File fastqc1_zip = "~{read1_name}_fastqc.zip"
    File fastqc2_html = "~{read2_name}_fastqc.html"
    File fastqc2_zip = "~{read2_name}_fastqc.zip"
    Int read1_seq = read_string("READ1_SEQS")
    Int read2_seq = read_string("READ2_SEQS")
    String read_pairs = read_string("READ_PAIRS")
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
task fastqc_se {
  input {
    File read1
    String read1_name = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    Int? cpus = 2
    String docker="quay.io/staphb/fastqc:0.11.9"
  }
  command <<<
    # capture date and version
    date | tee DATE
    fastqc --version | grep FastQC | tee VERSION

    fastqc --outdir $PWD --threads ~{cpus} ~{read1}

    unzip -p ~{read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS

    READ_SEQS=$(unzip -p ~{read1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
  >>>
  output {
    File fastqc_html = "~{read1_name}_fastqc.html"
    File fastqc_zip = "~{read1_name}_fastqc.zip"
    Int number_reads = read_string("READ1_SEQS")
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker:  "~{docker}"
    memory:  "4 GB"
    cpu:   2
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}
task quast {
  input {
    File assembly
    String samplename
    String docker="quay.io/staphb/quast:5.0.2"
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
    CODE

  >>>
  output {
    File quast_report = "${samplename}_report.tsv"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    Int genome_length = read_int("GENOME_LENGTH")
    Int number_contigs = read_int("NUMBER_CONTIGS")
    Int n50_value = read_int("N50_VALUE")
  }
  runtime {
    docker:  "~{docker}"
    memory:  "2 GB"
    cpu:   2
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}
task cg_pipeline {
  input {
    File  read1
    File?  read2
    String  samplename
    String  docker="quay.io/staphb/lyveset:1.1.4f"
    String  cg_pipe_opts="--fast"
    Int genome_length
  }
  command <<<
    # date and version control
    date | tee DATE

    run_assembly_readMetrics.pl ~{cg_pipe_opts} ~{read1} ~{read2} -e ~{genome_length} > ~{samplename}_readMetrics.tsv
    
    
    python3 <<CODE
    import csv
    #grab output average quality and coverage scores by column header
    with open("~{samplename}_readMetrics.tsv",'r') as tsv_file:
      tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
      for line in tsv_reader:
          if "_1" or "_R1" in line["File"]:
            with open("R1_MEAN_Q", 'wt') as r1_mean_q:
              r1_mean_q.write(line["avgQuality"])
            coverage = float(line["coverage"])
            print(coverage)
          if "_2" or "_R2" in line["File"]:
            with open("R2_MEAN_Q", 'wt') as r2_mean_q:
              r2_mean_q.write(line["avgQuality"])
            coverage += float(line["coverage"])
            print()
            with open("EST_COVERAGE", 'wt') as est_coverage:
              est_coverage.write(str(coverage))
    CODE

  >>>
  output {
    File  cg_pipeline_report = "${samplename}_readMetrics.tsv"
    String  cg_pipeline_docker   = docker
    String  pipeline_date = read_string("DATE")
    Float r1_mean_q = read_float("R1_MEAN_Q")
    Float? r2_mean_q = read_float("R2_MEAN_Q")
    Float est_coverage = read_float("EST_COVERAGE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
