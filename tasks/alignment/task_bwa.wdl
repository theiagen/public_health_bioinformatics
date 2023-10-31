version 1.0

task bwa {
  input {
    File read1
    File? read2
    String samplename
    File? reference_genome
    Int cpu = 6
    Int memory = 16
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
  }
  command <<<
    # date and version control
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version )" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # set reference genome
    if [[ ! -z "~{reference_genome}" ]]; then
      echo "User reference identified; ~{reference_genome} will be utilized for alignement"
      ref_genome="~{reference_genome}"
      bwa index "~{reference_genome}"
      # move to primer_schemes dir; bwa fails if reference file not in this location
    else
      ref_genome="/artic-ncov2019/primer_schemes/nCoV-2019/V3/nCoV-2019.reference.fasta"  
    fi

    # Map with BWA MEM; pipe to samtools sort to write sorted SAM file
    bwa mem \
      -t ~{cpu} \
      "${ref_genome}" \
      ~{read1} \
      ~{read2} | \
    samtools sort -@ ~{cpu} > ~{samplename}.sorted.sam
    
    # convert SAM to BAM that only includes aligned reads
    samtools view \
      -@ ~{cpu} \
      -F 4 \
      -b \
      -o ~{samplename}.sorted.bam \
      ~{samplename}.sorted.sam

    # convert SAM to BAM that only includes unaligned reads
    samtools view \
      -@ ~{cpu} \
      -f 4 \
      -b \
      -o ~{samplename}.sorted.unaligned-reads.bam \
      ~{samplename}.sorted.sam

    # if read2 was provided by user, extract both read1 and read2 from aligned and unaligned BAMs
    if [[ ! -z "~{read2}" ]]; then
      echo "Generating FASTQs for aligned and unaligned paired reads"
      samtools fastq \
        -F 4 \
        -1 ~{samplename}_R1.fastq.gz \
        -2 ~{samplename}_R2.fastq.gz \
        ~{samplename}.sorted.bam
      # note the lowercase 'f' here is imporant
      samtools fastq \
        -f 4 \
        -1 ~{samplename}_unaligned_R1.fastq.gz \
        -2 ~{samplename}_unaligned_R2.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    else
      echo "Generating FASTQs for aligned and unaligned single-end reads"
      samtools fastq \
        -@ ~{cpu} \
        -F 4 \
        -1 ~{samplename}_R1.fastq.gz \
        ~{samplename}.sorted.bam 
      # again, lowercase 'f' is important for getting all unaligned reads
      samtools fastq \
        -@ ~{cpu} \
        -f 4 \
        -1 ~{samplename}_unaligned_R1.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    fi

    # index BAMs
    samtools index ~{samplename}.sorted.bam 
    samtools index ~{samplename}.sorted.unaligned-reads.bam
  >>>
  output {
    String bwa_version = read_string("BWA_VERSION")
    String sam_version = read_string("SAMTOOLS_VERSION")
    File sorted_bam = "${samplename}.sorted.bam"
    File sorted_bai = "${samplename}.sorted.bam.bai"
    File read1_aligned = "~{samplename}_R1.fastq.gz"
    File? read2_aligned = "~{samplename}_R2.fastq.gz"
    File read1_unaligned = "~{samplename}_unaligned_R1.fastq.gz"
    File? read2_unaligned = "~{samplename}_unaligned_R2.fastq.gz"
    File sorted_bam_unaligned = "~{samplename}.sorted.unaligned-reads.bam"
    File sorted_bam_unaligned_bai = "~{samplename}.sorted.unaligned-reads.bam.bai"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    #maxRetries: 3
  }
}