version 1.0

task fastp_pe {
  input {
    File read1
    File read2
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int fastp_minlen = 75
    Int fastp_window_size=4
    Int fastp_quality_trim_score=30
    # adapter trimming and quality filtering are enabled by default, the flags below disable these functions to match trimmomatic
    String fastp_args = "--disable_adapter_trimming --disable_quality_filtering --disable_trim_poly_g --dont_eval_duplication"
    Int threads = 4
    Int disk_size = 100
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} --in2 ~{read2} \
    --out1 ~{samplename}_1P.fastq.gz --out2 ~{samplename}_2P.fastq.gz \
    --cut_right --cut_window_size ~{fastp_window_size} --cut_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_args} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task fastp_se {
  input {
    File read1
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int fastp_minlen = 75
    Int fastp_window_size=4
    Int fastp_quality_trim_score=30
    Int threads = 4
    # adapter trimming and quality filtering are enabled by default, the flags below disable these functions
    String fastp_params = "--disable_adapter_trimming --disable_quality_filtering"
    Int disk_size = 100
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} \
    --out1 ~{samplename}_1P.fastq.gz \
    --cut_right --cut_window_size ~{fastp_window_size} --cut_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_params} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task fastp_pe_theiaprok {
  input {
    File read1
    File read2
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int disk_size = 100
    Int fastp_window_size = 20
    Int fastp_quality_trim_score = 30
    Int fastp_minlen = 50
    # -g enables polyg trimming with default value of 10
    String fastp_args = "--detect_adapter_for_pe -g -5 20 -3 20"
    Int threads = 4
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} --in2 ~{read2} \
    --out1 ~{samplename}_1P.fastq.gz --out2 ~{samplename}_2P.fastq.gz \
    --unpaired1 ~{samplename}_1U.fastq.gz --unpaired2 ~{samplename}_2U.fastq.gz \
    --cut_right --cut_right_window_size ~{fastp_window_size} --cut_right_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_args} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File read2_trimmed = "~{samplename}_2P.fastq.gz"
    File read1_trimmed_unpaired = "~{samplename}_1U.fastq.gz"
    File read2_trimmed_unpaired = "~{samplename}_2U.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "quay.io/staphb/fastp:0.23.2"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}

task fastp_se_theiaprok {
  input {
    File read1
    String samplename
    String docker = "quay.io/staphb/fastp:0.23.2"
    Int disk_size = 100
    Int fastp_window_size = 20
    Int fastp_quality_trim_score = 30
    Int fastp_minlen = 50
    # -g enables polyg trimming with default value of 10
    # --detect_adapter_for_pe argument was removed 
    String fastp_args = "-g -5 20 -3 20"
    Int threads = 4
  }
  command <<<
    # date 
    date | tee DATE

    fastp \
    --in1 ~{read1} \
    --out1 ~{samplename}_1P.fastq.gz \
    --cut_right --cut_right_window_size ~{fastp_window_size} --cut_right_mean_quality ~{fastp_quality_trim_score} \
    --length_required ~{fastp_minlen} \
    --thread ~{threads} \
    ~{fastp_args} \
    --html ~{samplename}_fastp.html --json ~{samplename}_fastp.json
  >>>
  output {
    File read1_trimmed = "~{samplename}_1P.fastq.gz"
    File fastp_stats = "~{samplename}_fastp.html"
    String version = "~{docker}"
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "quay.io/staphb/fastp:0.23.2"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}