version 1.0

task trimmomatic_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/trimmomatic:0.40"
    Int trimmomatic_min_length = 75
    Int trimmomatic_window_size = 4
    Int trimmomatic_quality_trim_score = 30
    Int cpu = 4
    Int disk_size = 100
    Int memory = 8
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    fi
    
    trimmomatic PE \
    ~{trimmomatic_args} \
    ~{read1} ~{read2} \
    -baseout ~{samplename}.fastq.gz \
    SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_quality_trim_score} \
    MINLEN:~{trimmomatic_min_length} &> ~{samplename}.trim.stats.txt

  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    String trimmomatic_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task trimmomatic_se {
  input {
    File read1
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/trimmomatic:0.40"
    Int trimmomatic_min_length = 25
    Int trimmomatic_window_size = 4
    Int trimmomatic_quality_trim_score = 30
    Int cpu = 4
    Int disk_size = 100
    Int memory = 8
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    trimmomatic SE \
    ~{trimmomatic_args} \
    ~{read1} \
    ~{samplename}_trimmed.fastq.gz \
    SLIDINGWINDOW:~{trimmomatic_window_size}:~{trimmomatic_quality_trim_score} \
    MINLEN:~{trimmomatic_min_length} > ~{samplename}.trim.stats.txt
  >>>
  output {
    File read1_trimmed = "~{samplename}_trimmed.fastq.gz"
    File trimmomatic_stats = "~{samplename}.trim.stats.txt"
    String version = read_string("VERSION")
    String pipeline_date = read_string("DATE")
    String trimmomatic_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}