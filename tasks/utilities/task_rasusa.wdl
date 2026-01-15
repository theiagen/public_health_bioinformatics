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
    # RASUA Parameters
    #  --bases [STRING] Explicitly set the number of bases required e.g., 4.3kb, 7Tb, 9000, 4.1MB. If this option is given, --coverage and --genome-size are ignored
    #  --coverage [FLOAT] The desired coverage to sub-sample the reads to. If --bases is not provided, this option and --genome-size are required
    #  --genome_length [STRING] Genome size to calculate coverage with respect to. e.g., 4.3kb, 7Tb, 9000, 4.1MB
    #  --seed [INTERGER] Random seed to use
    #  --frac [FLOAT] Subsample to a fraction of the reads - e.g., 0.5 samples half the reads
    #  --num [INTEGER] Subsample to a specific number of reads
    String? bases
    Float coverage = 250
    String genome_length
    Int? seed
    Float? frac
    Int? num 
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
    # ignore coverage and genome length if frac input provided
    if [ -z "~{frac}" ]; then
      COVERAGE="--coverage ~{coverage} --genome-size ~{genome_length}"
    else
      COVERAGE=""
    fi

    # run rasusa for read sampling
    rasusa reads \
      ${COVERAGE} \
      ~{'--seed ' + seed} \
      ~{'--bases ' + bases} \
      ~{'--frac ' + frac} \
      ~{'--num ' + num} \
      ${OUTPUT_FILES} \
      ~{read1} ~{read2}
  >>>
  output {
    File read1_subsampled = "~{samplename}_subsampled_R1.fastq.gz"
    File? read2_subsampled = "~{samplename}_subsampled_R2.fastq.gz"
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
