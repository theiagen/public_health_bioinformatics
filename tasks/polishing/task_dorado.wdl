version 1.0

task dorado {
  input {
    File read1
    File unpolished_fasta
    String samplename

    String model

    Int cpu
    Int memory
    Int disk_size
    String docker
}
command <<<
  set -euo pipefail
  # align with dorado aligner before polishing

  dorado aligner \
    -t ~{cpu} \
    ~{unpolished_fasta} \
    ~{read1} > ~{samplename}.bam

  samtools sort \
    -@ ~{cpu} \
    ~{samplename}.bam > ~{samplename}_sorted.bam
  samtools index ~{samplename}_sorted.bam

  # figure out how to sort out basecalling model identification
  #
  dorado polish \
     -t ~{cpu} \
     --models-directory \
     ~{samplename}_sorted.bam \
     ~{unpolished_fasta} > ~{samplename}_polished.fasta


>>>
