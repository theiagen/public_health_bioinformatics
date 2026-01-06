version 1.0

task bbduk {
  input {
    File read1
    File read2
    String samplename
    Int memory = 8
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bbtools:38.76"
    Int disk_size = 100
    File? adapters_fasta
    File? phix_fasta
  }
  command <<<
    # date and version control
    date | tee DATE

    # Repairing disordered reads (if they exist) so that the first read in file 1 is the same mate of the first read in file 2
    echo "Repairing paired-end reads to ensure correct order..."
    repair.sh \
      in=~{read1} \
      in2=~{read2} \
      out=~{samplename}.raw_1.fastq.gz \
      out2=~{samplename}.raw_2.fastq.gz

    # Set phix fasta
    if [[ -n "~{phix_fasta}" ]]; then
      phix_fasta="~{phix_fasta}"
      echo "Using user supplied FASTA file for phiX: '~{phix_fasta}'"
    else
      phix_fasta="/bbmap/resources/phix174_ill.ref.fa.gz"
      echo "Using default phiX FASTA file: '/bbmap/resources/phix174_ill.ref.fa.gz'"
    fi

    # Attempt phiX removal first to remove contamination
    echo "Filtering and removing reads contaminated with phiX..."
    bbduk.sh \
      in=~{samplename}.raw_1.fastq.gz \
      in2=~{samplename}.raw_2.fastq.gz \
      out=~{samplename}.rm_phix_1.fastq.gz \
      out2=~{samplename}.rm_phix_2.fastq.gz \
      ref=${phix_fasta} \
      stats=~{samplename}.phix.stats.txt \
      statscolumns=5 \
      k=31 hdist=1 ordered=t

    # Set adapter fasta
    if [[ -n "~{adapters_fasta}" ]]; then
      adapter_fasta="~{adapters_fasta}"
      echo "Using user supplied FASTA file for adapters: '~{adapters_fasta}'"
    else
      adapter_fasta="/bbmap/resources/adapters.fa"
      echo "Using default adapters FASTA file: '/bbmap/resources/adapters.fa'"
    fi

    # Trim adapters
    echo "Trimming adapters and capturing those reads..."
    bbduk.sh \
      in=~{samplename}.rm_phix_1.fastq.gz \
      in2=~{samplename}.rm_phix_2.fastq.gz \
      out=~{samplename}.rm_adpt_1.fastq.gz \
      out2=~{samplename}.rm_adpt_2.fastq.gz \
      ref=${adapter_fasta} \
      stats=~{samplename}.adapters.stats.txt \
      statscolumns=5 \
      k=23 ktrim=r mink=11 hdist=1 tpe=t tbo=t ordered=t

    # Rename output files to final cleaned filenames
    mv ~{samplename}.rm_adpt_1.fastq.gz ~{samplename}_1.clean.fastq.gz
    mv ~{samplename}.rm_adpt_2.fastq.gz ~{samplename}_2.clean.fastq.gz

  >>>
  output {
    File read1_clean = "${samplename}_1.clean.fastq.gz"
    File read2_clean = "${samplename}_2.clean.fastq.gz"
    File adapter_stats = "${samplename}.adapters.stats.txt"
    File phiX_stats = "${samplename}.phix.stats.txt"
    String bbduk_docker = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}

task bbduk_se {
  input {
    File read1_trimmed
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bbtools:38.76"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
    File? adapters_fasta
    File? phix_fasta
  }
  command <<<
    # date and version control
    date | tee DATE

    # set adapter fasta
    if [[ ! -z "~{adapters_fasta}" ]]; then
      echo "Using user supplied FASTA file for adapters..."
      adapter_fasta="~{adapters_fasta}"
    else
      echo "User did not supply adapters FASTA file, using default adapters.fa file..."
      adapter_fasta="/bbmap/resources/adapters.fa" 
    fi

    # set phix fasta
    if [[ ! -z "~{phix_fasta}" ]]; then
      echo "Using user supplied FASTA file for phiX..."
      phix_fasta="~{phix_fasta}"
    else
      echo "User did not supply phiX FASTA file, using default phix174_ill.ref.fa.gz file..."
      phix_fasta="/bbmap/resources/phix174_ill.ref.fa.gz"
    fi

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}.rmadpt_1.fastq.gz ref=${adapter_fasta} stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo ordered=t

    bbduk.sh in1=~{read1_trimmed} out1=~{samplename}_1.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=${phix_fasta} k=31 hdist=1 stats=~{samplename}.phix.stats.txt ordered=t
  >>>
  output {
    File read1_clean = "${samplename}_1.clean.fastq.gz"
    File adapter_stats = "${samplename}.adapters.stats.txt"
    File phiX_stats = "${samplename}.phix.stats.txt"
    String bbduk_docker   = docker
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}