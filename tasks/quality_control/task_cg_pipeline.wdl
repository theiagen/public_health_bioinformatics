version 1.0

task cg_pipeline {
  input {
    File read1
    File? read2
    String samplename
    String docker="quay.io/staphb/lyveset:1.1.4f"
    String cg_pipe_opts="--fast"
    Int genome_length
  }
  command <<<
    # date and version control
    date | tee DATE

    run_assembly_readMetrics.pl ~{cg_pipe_opts} ~{read1} ~{read2} -e ~{genome_length} > ~{samplename}_readMetrics.tsv
    
    python3 <<CODE
    import csv
    #grab output average quality and coverage scores by column header
    coverage = 0.0
    with open("~{samplename}_readMetrics.tsv",'r') as tsv_file:
      tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
      for line in tsv_reader:
        if "~{read1}" in line["File"]:
          with open("R1_MEAN_Q", 'wt') as r1_mean_q:
            r1_mean_q.write(line["avgQuality"])
          
          # run_assembly_readMetrics can report coverage as '.'
          try:
            coverage = float(line["coverage"])
          except ValueError:
            continue
          print(coverage)
          
        else:
          with open("R2_MEAN_Q", 'wt') as r2_mean_q:
            r2_mean_q.write(line["avgQuality"])
          # run_assembly_readMetrics can report coverage as '.'
          try:
            coverage += float(line["coverage"])
          except ValueError:
            continue

      with open("EST_COVERAGE", 'wt') as est_coverage:
        est_coverage.write(str(coverage))
    CODE

    # R2_MEAN_Q to make SE workflow work otherwise read_float fails
    if [[ ! -f R2_MEAN_Q ]] ; then
      echo "0.0" > R2_MEAN_Q
    fi

  >>>
  output {
    File cg_pipeline_report = "${samplename}_readMetrics.tsv"
    String cg_pipeline_docker   = docker
    String pipeline_date = read_string("DATE")
    Float r1_mean_q = read_float("R1_MEAN_Q")
    Float r2_mean_q = read_float("R2_MEAN_Q")
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