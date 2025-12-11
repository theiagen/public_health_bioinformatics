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
    Int? trimmomatic_base_crop
    Int cpu = 4
    String trimmomatic_args = "-phred33"
    Int disk_size = 100
    Int memory = 8
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    CROPPING_VAR=""
    # if trimmomatic base chop is defined (-n means not empty), determine average readlength of the input reads
    if [ -n "~{trimmomatic_base_crop}" ]; then
      # determine the average read length of the input reads
      read_length_r1=$(zcat ~{read1} | awk '{if(NR%4==2) {bases+=length($0)} } END {print bases/(NR/4)}')
      read_length_r2=$(zcat ~{read2} | awk '{if(NR%4==2) {bases+=length($0)} } END {print bases/(NR/4)}')

      # take the average of the two read lengths without using bc and remove the end base chop
      avg_readlength=$(python3 -c "print(int(((${read_length_r1} + ${read_length_r2}) / 2) - ~{trimmomatic_base_crop}))")
    
      # HEADCROP: number of bases to remove from the start of the read
      # CROP: number of bases to KEEP, from the start of the read
      CROPPING_VAR="HEADCROP:~{trimmomatic_base_crop} CROP:$avg_readlength"
    fi
    
    trimmomatic PE \
    ~{trimmomatic_args} \
    -threads ~{cpu} \
    ~{read1} ~{read2} \
    -baseout ~{samplename}.fastq.gz \
    "${CROPPING_VAR}" \
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
    String trimmomatic_args = "-phred33"
    Int disk_size = 100
    Int memory = 8
  }
  command <<<
    # date and version control
    date | tee DATE
    trimmomatic -version > VERSION && sed -i -e 's/^/Trimmomatic /' VERSION

    trimmomatic SE \
    ~{trimmomatic_args} \
    -threads ~{cpu} \
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