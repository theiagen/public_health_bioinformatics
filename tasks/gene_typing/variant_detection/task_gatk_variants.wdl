version 1.0

task gatk_variants {
  input {
    File bam
    File bai
    File reference_genome
    File? intervals_file
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/gatk:4.6.2.0-dev"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100


    Int map_qual = 60 # set to mirror v4.6.0 default
    Int base_quality = 13 # set to mirror v4.6.0 default
    Int min_coverage = 10 # set to mirror v4.6.0 default
    Float min_frac = 0 # set to mirror v4.6.0 default
    Int min_quality = 100 # set to mirror v4.6.0 default
    Int maxsoft = 10 # set to mirror v4.6.0 default
  }
  command <<<
    # fail hard
    set -euo pipefail

    # obtain version
    gatk --version | grep GATK | grep -Po "v[^ ]+$" | tee VERSION

    # call HaplotypeCaller
    gatk \
      --java-options "-Xms~{memory}G -Xmx~{memory}G" \
      HaplotypeCaller \
      -R ~{reference_genome} \
      -I ~{bam} \
      -o ~{samplename}_haplotypecall.g.vcf.gz
      ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
      --tmp-dir .

    # call GenotypeGVCFs
    gatk --java-options "-Xmx~{memory}G" \
      GenotypeGVCFs \
      -R ~{reference_genome} \
      -V ~{samplename}_haplotypecall.g.vcf.gz \
      ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
      -O ~{samplename}_genotype.g.vcf.gz \

    # call VariantFiltration
    gatk --java-options "-Xmx~{memory}G" VariantFiltration \
      -R ~{reference_genome} \
      -V ~{samplename}_genotype.g.vcf.gz \
      -O ~{samplename}_filtered.g.vcf.gz

    # call SelectVariants
    gatk --java-options "-Xmx~{memory}G" SelectVariants \
      -V ~{samplename}_filtered.g.vcf.gz \
      -O ~{samplename}_selected.g.vcf.gz
  >>>
  output {
    String gatk_version = read_string("VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: cpu
      disks: "local-disk " + disk_size + " SSD"
      preemptible: 0
      maxRetries: 3
  }
}
