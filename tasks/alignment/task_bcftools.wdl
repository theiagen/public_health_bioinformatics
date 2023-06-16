version 1.0

task bcftools_mpileup {
  meta {
    description: "Build pileup file with bcftools"
  }
  input {
    File bam
    File bai
    File reference
    String samplename
    String docker = "staphb/bcftools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    bcftools mpileup -Ou -f ~{reference} ~{bam} | bcftools call -Ou -mv | bcftools norm -f ~{reference} -Oz -o "~{samplename}".vcf.gz
    bcftools query -f'%CHROM\t%POS0\t%END\n' "~{samplename}".vcf.gz > "~{samplename}".variants.bed
  >>>
  output {
    File bcftools_vcf = "~{samplename}.vcf.gz"
    File bcftools_variants_bed = "~{samplename}.variants.bed"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task bedtools_mask {
  meta {
    description: "Build bedfile with sites to be masked with bedtools"
  }
  input {
    File bam
    File bai
    File variants
    String samplename
    String docker = "staphb/bedtools:2.30.0"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    bedtools genomecov -bga -ibam ~{bam} | awk '$4 < 10' > low_coverage_sites.bed
    bedtools subtract -a low_coverage_sites.bed -b ~{variants} > "~{samplename}".mask.bed
  >>>
  output {
    File bedtools_mask_bed = "~{samplename}.mask.bed"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}

task bcftools_consensus {
  meta {
    description: "Build consensus sequence file with bcftools"
  }
  input {
    File vcf
    File reference
    File mask
    String samplename
    String docker = "staphb/bcftools:1.17"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    bcftools index -t ~{vcf}
    bcftools consensus -f ~{reference} --mask ~{mask} --mask-with N --mark-del - ~{vcf} > "~{samplename}"_consensus.fasta
    sed -i -e 's/^>\(.*\)/>~{samplename}/' "~{samplename}"_consensus.fasta
    # remove non-ATCGN bases and replace them with N
    sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D|-]/N/g' "~{samplename}"_consensus.fasta > "~{samplename}"_consensus.fasta

  >>>
  output {
    File bcftools_consensus_fasta = "~{samplename}_consensus.fasta"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}