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
    set -euo pipefail
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

    # set cat command based on compression
    if [[ "~{read1}" == *".gz" ]] ; then
      cat_reads="zcat"
    else
      cat_reads="cat"
    fi

    echo -e "\ninput R1 has $(${cat_reads} ~{read1} | grep -c '^@') reads as input"
    echo "input R2 has $(${cat_reads} ~{read2} | grep -c '^@') reads as input"

    # Map with BWA MEM; pipe to samtools sort to write sorted SAM file
    bwa mem \
      -t ~{cpu} \
      "${ref_genome}" \
      ~{read1} \
      ~{read2} | \
    samtools sort \
      -@ ~{cpu} - \
      > ~{samplename}.sorted.sam
    
    # convert SAM to BAM that does not include unaligned reads; "-F 4" = exclude unaligned reads
    # FYI - secondary and supplementary alignments are included in output BAM file
    samtools view \
      -@ ~{cpu} \
      -F 4 \
      -b \
      -o ~{samplename}.sorted.aligned.bam \
      ~{samplename}.sorted.sam

    # create a separate BAM that only includes aligned reads (no secondary or supplementary alignments) for creating aligned FASTQ files
    ## 0x904 = 4 (unaligned) + 100 (secondary) + 800 (supplementary)
    # FYI the default filtering is 0x900 (secondary and supplemental) so when we previously specified -F 4, we were only filtering out unaligned reads, and keeping secondary and supplemental alignments
    # this update to 0x904 means we are now filtering out unaligned, secondary, and supplemental alignments
    samtools view \
      -@ ~{cpu} \
      -F 0x904 \
      -b \
      -o ~{samplename}.sorted.aligned-filtered.bam \
      ~{samplename}.sorted.sam

    # convert SAM to BAM that only includes unaligned reads
    samtools view \
      -@ ~{cpu} \
      -f 4 \
      -b \
      -o ~{samplename}.sorted.unaligned-reads.bam \
      ~{samplename}.sorted.sam

    # see here for "samtools fastq" options: https://www.htslib.org/doc/samtools-fasta.html
    # TL;DR is that "samtools fastq -1 R1.fastq -2 R2.fastq" works with paired-end inputs and will output R1 and R2 reads to separate files due to tags in the SAM & BAM file

    # AFAIK for single end alignments w/ bwa mem, the output SAM/BAM files do not have tags to differentiate between R1 and R2 reads, so the "samtools fastq -0 R1.fastq" command is used to output all reads to a single file

    # if read2 was provided by user, extract both read1 and read2 from aligned and unaligned BAMs
    if [[ ! -z "~{read2}" ]]; then
      echo -e "\nGenerating FASTQs for aligned reads"
      samtools fastq \
        -@ ~{cpu} \
        -F 4 \
        -1 ~{samplename}_R1.fastq.gz \
        -2 ~{samplename}_R2.fastq.gz \
        ~{samplename}.sorted.aligned-filtered.bam
      echo "Generating FASTQs for unaligned reads"
      # note the lowercase 'f' here is imporant
      samtools fastq \
        -@ ~{cpu} \
        -f 4 \
        -1 ~{samplename}_unaligned_R1.fastq.gz \
        -2 ~{samplename}_unaligned_R2.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    else
      echo -e "\nGenerating FASTQs for aligned single-end reads\n"
      samtools fastq \
        -@ ~{cpu} \
        -F 4 \
        -0 ~{samplename}_R1.fastq.gz \
        ~{samplename}.sorted.aligned-filtered.bam
      echo -e "Generating FASTQs for unaligned single-end reads\n"
      # again, lowercase 'f' is important for getting all unaligned reads
      samtools fastq \
        -@ ~{cpu} \
        -f 4 \
        -0 ~{samplename}_unaligned_R1.fastq.gz \
        ~{samplename}.sorted.unaligned-reads.bam
    fi

    # index BAMs
    samtools index ~{samplename}.sorted.aligned.bam 
    samtools index ~{samplename}.sorted.unaligned-reads.bam

    # count output reads to ensure we are outputting all reads, regardless if the aligned or not
    # if read2 does exist as input, count both R1 and R2
    if [[ ! -z "~{read2}" ]]; then
      echo -e "\noutput R1_aligned has $(zcat ~{samplename}_R1.fastq.gz | grep -c '^@') reads as input"
      echo "output R2_aligned has $(zcat ~{samplename}_R2.fastq.gz | grep -c '^@') reads as input"
      echo
      echo "output R1_unaligned has $(zcat ~{samplename}_unaligned_R1.fastq.gz | grep -c '^@') reads as input"
      echo "output R2_unaligned has $(zcat ~{samplename}_unaligned_R2.fastq.gz | grep -c '^@') reads as input"
    # else = if read2 does not exist as input, only count R1
    else
      echo "output R1_aligned has $(zcat ~{samplename}_R1.fastq.gz | grep -c '^@') reads as input"
      echo
      echo "output R1_unaligned has $(zcat ~{samplename}_unaligned_R1.fastq.gz | grep -c '^@') reads as input"
    fi
  >>>
  output {
    String bwa_version = read_string("BWA_VERSION")
    String sam_version = read_string("SAMTOOLS_VERSION")
    File sorted_bam = "${samplename}.sorted.aligned.bam"
    File sorted_bai = "${samplename}.sorted.aligned.bam.bai"
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
    maxRetries: 3
  }
}

task bwa_all {
  input {
    File draft_assembly_fasta
    File read1
    File read2
    String samplename

    Int cpu = 6
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bwa:0.7.18"
    Int memory = 16
  }
  command <<<
    set -euo pipefail

    # Get version
    echo "BWA $(bwa 2>&1 | grep Version )" | tee BWA_VERSION

    if [[ ! -f "~{draft_assembly_fasta}.bwt" ]]; then
      echo "Indexing reference genome: ~{draft_assembly_fasta}"
      bwa index ~{draft_assembly_fasta}
    else
      echo "Reference genome is already indexed: ~{draft_assembly_fasta}"
    fi
    
    bwa mem -t ~{cpu} -a ~{draft_assembly_fasta} ~{read1} > ~{samplename}_R1.sam
    bwa mem -t ~{cpu} -a ~{draft_assembly_fasta} ~{read2} > ~{samplename}_R2.sam

  >>>
  output {
    File read1_sam = "~{samplename}_R1.sam"
    File read2_sam = "~{samplename}_R2.sam"
    String bwa_version = read_string("BWA_VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
