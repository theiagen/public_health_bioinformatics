version 1.0

task rasusa {
  meta {
    description: "Randomly subsample sequencing reads to a specified coverage (https://github.com/mbhall88/rasusa)"
    volatile: true
  }
  input {
    File read1
    File? read2
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/rasusa:2.1.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 8
    # RASUSA Parameters
    #  --bases [STRING] Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB. If this option is given, --coverage and --genome-size are ignored
    #  --coverage [FLOAT] The desired coverage to sub-sample the reads to. If --bases is not provided, this option and --genome-size are required
    #  --genome_length [STRING] Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB
    #  --seed [INTERGER] Random seed to use
    #  --frac [FLOAT] Subsample to a fraction of the reads - e.g., 0.5 samples half the reads
    #  --num [INTEGER] Subsample to a specific number of reads
    String? num_bases
    Float coverage = 250
    String? genome_length
    Int? seed
    Float? fraction_of_reads
    Int? num_reads
  }
  command <<<
    # fail hard
    set -euo pipefail
    rasusa --version | tee VERSION
    # set single-end or paired-end outputs
    if [ -z "~{read2}" ]; then
      OUTPUT_FILES="-o ~{samplename}_subsampled_R1.fastq.gz"
    else
      OUTPUT_FILES="-o ~{samplename}_subsampled_R1.fastq.gz -o ~{samplename}_subsampled_R2.fastq.gz"
    fi

    # valid options: (coverage + genome_length), fraction_of_reads, num_bases, or num_reads
    # fraction_of_reads, num_bases, and num_reads should be mutually exclusive and take precedence over coverage/genome_length
    OVERRIDE_PARAMS=0
    if [ -n "~{fraction_of_reads}" ]; then
      OVERRIDE_PARAMS=$((OVERRIDE_PARAMS + 1));
    fi
    if [ -n "~{num_bases}" ]; then
      OVERRIDE_PARAMS=$((OVERRIDE_PARAMS + 1));
    fi
    if [ -n "~{num_reads}" ]; then
      OVERRIDE_PARAMS=$((OVERRIDE_PARAMS + 1));
    fi

    if [ "$OVERRIDE_PARAMS" -gt 1 ]; then
      echo "ERROR: fraction_of_reads, num_bases, and num_reads are mutually exclusive; specify only one"
      exit 1
    elif [ "$OVERRIDE_PARAMS" -eq 1 ]; then
      echo "INFO: Using fraction_of_reads, num_bases, or num_reads parameter; ignoring coverage and genome_length"
      COVERAGE=""
    elif [ -n "~{genome_length}" ] && [ -n "~{coverage}" ]; then
      COVERAGE="--coverage ~{coverage} --genome-size ~{genome_length}"
    else
      echo "ERROR: Provide genome_length and coverage, or one of: fraction_of_reads, num_bases, num_reads"
      exit 1
    fi

    # run rasusa for read sampling
    rasusa reads -v \
      ${COVERAGE} \
      ~{'--seed ' + seed} \
      ~{'--bases ' + num_bases} \
      ~{'--frac ' + fraction_of_reads} \
      ~{'--num ' + num_reads} \
      ${OUTPUT_FILES} \
      ~{read1} ~{read2} \
      2> rasusa.log
  >>>
  output {
    File read1_subsampled = "~{samplename}_subsampled_R1.fastq.gz"
    File? read2_subsampled = "~{samplename}_subsampled_R2.fastq.gz"
    File rasusa_log = "rasusa.log"
    String rasusa_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
