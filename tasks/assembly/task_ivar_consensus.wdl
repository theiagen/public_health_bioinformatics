version 1.0

task consensus {
  input {
    File bamfile
    String samplename
    File? reference_genome
    Boolean count_orphans = true
    Int max_depth = "600000"
    Boolean disable_baq = true
    Boolean all_positions = false
    Int min_bq = "0"
    Int min_qual = "20"
    Float? consensus_min_freq 
    Int? consensus_min_depth
    String char_unknown = "N"
    Boolean skip_N = false
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
  }
  command <<<
    #set -euo pipefail to avoid silent failure
    set -euo pipefail
    # date and version control
    date | tee DATE
    ivar version | head -n1 | tee IVAR_VERSION
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

    # index bam file
    samtools index ~{bamfile}

    # create fastas for each contig in reference genome file
    sed -E 's/^>([^ ]+).*$/>\1/' ~{reference_genome} | \
      awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".refcontig.fasta"} {print >> F}'

    for contig in *.refcontig.fasta; do
      contig_name=$(basename $contig ".refcontig.fasta")
      echo "DEBUG: Extracting alignments to reference contig: $contig_name"
      samtools view ~{bamfile} ${contig_name} -b -o ${contig_name}.bam

      if [[ ! -s ${contig_name}.bam ]]; then
        echo "WARNING: No reads aligned to reference contig ${contig_name}. Skipping consensus generation for this contig."
        continue
      fi

      # call mpileup
      samtools mpileup \
        ~{true = "--count-orphans" false = "" count_orphans} \
        -d ~{max_depth} \
        ~{true = "--no-BAQ" false = "" disable_baq} \
        -Q ~{min_bq} \
        --reference ${contig} \
        ~{true = "-aa" false = "" all_positions} \
        ${contig_name}.bam \
        > ${contig_name}.mpileup

      # call consensus
      cat ${contig_name}.mpileup | \
      ivar consensus \
        -p ${contig_name}.consensus \
        -q ~{min_qual} \
        -t ~{consensus_min_freq} \
        -m ~{consensus_min_depth} \
        -n ~{char_unknown} \
        ~{true = "-k" false = "" skip_N}

      # clean up fasta header
      echo ">${contig_name}" > ${contig_name}.consensus.fasta
      grep -v ">" ${contig_name}.consensus.fa >> ${contig_name}.consensus.fasta
    done

    # Combine all consensus sequences into a single fasta file
    cat *.consensus.fasta > ~{samplename}.ivar.consensus.fasta
    cat *.mpileup > ~{samplename}.mpileup
  >>>
  output {
    File consensus_seq = "~{samplename}.ivar.consensus.fasta"
    File sample_mpileup = "~{samplename}.mpileup"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String pipeline_date = read_string("DATE")
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
