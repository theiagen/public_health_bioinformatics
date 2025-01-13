version 1.0

task cg_pipeline {
  input {
    File read1
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/lyveset:1.1.4f"
    String cg_pipe_opts = "--fast"
    Int genome_length
    Int disk_size = 50
    Int memory = 2 # added a default value here
    Int cpu = 1 # added a default value here
  }
  command <<<
    # date and version control
    date | tee DATE

    run_assembly_readMetrics.pl ~{cg_pipe_opts} ~{read1} ~{read2} -e ~{genome_length} > ~{samplename}_readMetrics.tsv

    # repeat for concatenated read file
    # run_assembly_readMetrics.pl extension awareness    
    if [[ "~{read1}" == *".gz" ]] ; then
      extension=".gz"
    else
      extension=""
    fi
    cat ~{read1} ~{read2} > ~{samplename}_concat.fastq"${extension}"
    run_assembly_readMetrics.pl ~{cg_pipe_opts} ~{samplename}_concat.fastq"${extension}" -e ~{genome_length} > ~{samplename}_concat_readMetrics.tsv
        
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
          with open("R1_MEAN_LENGTH", 'wt') as r1_mean_length:
            r1_mean_length.write(line["avgReadLength"])       

          # run_assembly_readMetrics can report coverage as '.'
          try:
            coverage = float(line["coverage"])
          except ValueError:
            continue
          print(coverage)
          
        else:
          with open("R2_MEAN_Q", 'wt') as r2_mean_q:
            r2_mean_q.write(line["avgQuality"])
          with open("R2_MEAN_LENGTH", 'wt') as r2_mean_length:
            r2_mean_length.write(line["avgReadLength"]) 
          # run_assembly_readMetrics can report coverage as '.'
          try:
            coverage += float(line["coverage"])
          except ValueError:
            continue

      with open("EST_COVERAGE", 'wt') as est_coverage:
        est_coverage.write(str(coverage))

    # parse concatenated read metrics
    # grab output average quality and coverage scores by column header
    with open("~{samplename}_concat_readMetrics.tsv",'r') as tsv_file_concat:
      tsv_reader_concat = list(csv.DictReader(tsv_file_concat, delimiter="\t"))
      for line in tsv_reader_concat:
        if "~{samplename}_concat" in line["File"]:
          with open("COMBINED_MEAN_Q", 'wt') as combined_mean_q:
            combined_mean_q.write(line["avgQuality"])
          with open("COMBINED_MEAN_LENGTH", 'wt') as combined_mean_length:
            combined_mean_length.write(line["avgReadLength"])            

    CODE

    # R2_MEAN_Q to make SE workflow work otherwise read_float fails
    if [[ ! -f R2_MEAN_Q ]] ; then
      echo "0.0" > R2_MEAN_Q
    fi
    # same for R2_MEAN_LENGTH
    if [[ ! -f R2_MEAN_LENGTH ]] ; then
      echo "0.0" > R2_MEAN_LENGTH
    fi    
  >>>
  output {
    File cg_pipeline_report = "${samplename}_readMetrics.tsv"
    String cg_pipeline_docker = docker
    String pipeline_date = read_string("DATE")
    Float r1_mean_q = read_float("R1_MEAN_Q")
    Float r2_mean_q = read_float("R2_MEAN_Q")
    Float combined_mean_q = read_float("COMBINED_MEAN_Q")
    Float r1_mean_readlength = read_float("R1_MEAN_LENGTH")
    Float r2_mean_readlength = read_float("R2_MEAN_LENGTH")
    Float combined_mean_readlength = read_float("COMBINED_MEAN_LENGTH")
    Float est_coverage = read_float("EST_COVERAGE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}