version 1.0

task bbduk {
  input {
    File read1
    File read2
    String samplename
    Int memory = 8
    Int cpu = 4
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/bbtools:39.38_python"
    Int disk_size = 100

    File? adapters_fasta
    File? phix_fasta
    File? primers_fasta
    String? primers_literal

    Int primers_restrict_trim_length = 5
    Int primers_hamming_distance = 1
    Boolean primers_mask_middle = false
    Boolean primers_reverse_complement = true
  }
  command <<<
    # date and version control
    date | tee DATE

    # set adapter fasta
    if [[ ! -z "~{adapters}" ]]; then
      echo "Using user supplied FASTA file for adapters..."
      adapter_fasta="~{adapters}"
    # Repairing disordered reads (if they exist) so that the first read in file 1 is the same mate of the first read in file 2
    echo "Repairing paired-end reads to ensure correct order..."
    repair.sh \
      in=~{read1} \
      in2=~{read2} \
      out=~{samplename}.raw_1.fastq.gz \
      out2=~{samplename}.raw_2.fastq.gz

    else
      echo "User did not supply adapters FASTA file, using default adapters.fa file..."
      adapter_fasta="/bbmap/resources/adapters.fa" 
    fi

    # set phix fasta
    if [[ ! -z "~{phix}" ]]; then
      echo "Using user supplied FASTA file for phiX..."
      phix_fasta="~{phix}"
    else
      echo "User did not supply phiX FASTA file, using default phix174_ill.ref.fa.gz file..."
      phix_fasta="/bbmap/resources/phix174_ill.ref.fa.gz"
    fi


    bbduk.sh in1=~{samplename}.paired_1.fastq.gz in2=~{samplename}.paired_2.fastq.gz out1=~{samplename}.rmadpt_1.fastq.gz out2=~{samplename}.rmadpt_2.fastq.gz ref=${adapter_fasta} stats=~{samplename}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo ordered=t

    bbduk.sh in1=~{samplename}.rmadpt_1.fastq.gz in2=~{samplename}.rmadpt_2.fastq.gz out1=~{samplename}_1.clean.fastq.gz out2=~{samplename}_2.clean.fastq.gz outm=~{samplename}.matched_phix.fq ref=${phix_fasta} k=31 hdist=1 stats=~{samplename}.phix.stats.txt ordered=t
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
    File? adapters
    File? phix
  }
  command <<<
    # date and version control
    date | tee DATE

    # set adapter fasta
    if [[ ! -z "~{adapters}" ]]; then
      echo "Using user supplied FASTA file for adapters..."
      adapter_fasta="~{adapters}"
    else
      echo "User did not supply adapters FASTA file, using default adapters.fa file..."
      adapter_fasta="/bbmap/resources/adapters.fa" 
    fi

    # set phix fasta
    if [[ ! -z "~{phix}" ]]; then
      echo "Using user supplied FASTA file for phiX..."
      phix_fasta="~{phix}"
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