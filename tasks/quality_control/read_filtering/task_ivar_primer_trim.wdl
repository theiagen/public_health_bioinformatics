version 1.0

task primer_trim {
  input {
    File bamfile
    String samplename
    File primer_bed
    Boolean keep_noprimer_reads = true
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
  }
  String primer_name = basename(primer_bed)
  command <<<
    set -euo pipefail

    # date and version control
    echo "~{primer_name}" | tee PRIMER_NAME
    date | tee DATE

    ivar version | head -n1 | tee IVAR_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # trimming primers
    ivar trim \
      ~{true="-e" false="" keep_noprimer_reads} \
      -i ~{bamfile} \
      -b ~{primer_bed} \
      -p ~{samplename}.primertrim | tee ~{samplename}.primertrim.log

    # sorting and indexing the trimmed bams
    samtools sort \
      ~{samplename}.primertrim.bam \
      -o ~{samplename}.primertrim.sorted.bam

    samtools index ~{samplename}.primertrim.sorted.bam

    # soft-fail signal for downstream tasks when no reads aligned
    # we aren't filtering secondary and supplementary alignments because those are included
    # downstream dependencies that convert to dependencies will have to account for this
    mapped_count=$(samtools view -c ~{samplename}.primertrim.sorted.bam)
    echo "DEBUG: ${mapped_count} mapped records"
    if [[ "${mapped_count}" -eq 0 ]]; then
      echo "ERROR: no alignments remain following trimming" | tee PRIMER_TRIM_WARNING
    else
      touch PRIMER_TRIM_WARNING
    fi

    PCT=$(grep "Trimmed primers from" ~{samplename}.primertrim.log | perl -lape 's/Trimmed primers from (\S+)%.*/$1/')
    echo $PCT
    if [[ $PCT = -* ]]; then echo 0; else echo $PCT; fi > IVAR_TRIM_PCT
  >>>
  output {
    File trimmed_bam = "~{samplename}.primertrim.bam"
    File trim_sorted_bam = "~{samplename}.primertrim.sorted.bam"
    File trim_sorted_bai = "~{samplename}.primertrim.sorted.bam.bai"
    File trim_log = "~{samplename}.primertrim.log"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String pipeline_date = read_string("DATE")
    Float primer_trimmed_read_percent = read_float("IVAR_TRIM_PCT")
    String primer_bed_name = read_string("PRIMER_NAME")
    String primer_trim_warning = read_string("PRIMER_TRIM_WARNING")
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
