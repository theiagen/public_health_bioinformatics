version 1.0

task gatk_variants {
  input {
    String samplename
    File bam
    File bai
    File reference_genome
    Int ploidy = 1 # integer indicating ploidy (N); default to haploid
    File? intervals_file

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/gatk:4.6.2.0"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # fail hard
    set -euo pipefail

    # obtain version
    gatk --version | grep GATK | grep -Po "v[^ ]+$" | tee VERSION

    # create reference index and dictionary (have to symlink in case the source is not writeable)
    local_ref=$(basename ~{reference_genome})
    ln -s ~{reference_genome} ${local_ref}
    # index reference FASTA
    samtools faidx ${local_ref}
    gatk CreateSequenceDictionary -R ${local_ref}

    # identify read groups in BAM
    gatk GetSampleName \
      -I ~{bam} \
      -O read_groups.txt || echo "WARNING: No read groups found in BAM, assuming a single read group"

    # HaploTypeCaller needs read groups, so assume a single read classification if they don't exist
    # Currently, multiple read groups will function appropriately if inputted directly via ~{bam}
    if [ ! -f read_groups.txt ]; then
      picard -Xmx~{memory}G AddOrReplaceReadGroups \
        INPUT=~{bam} \
        OUTPUT=~{samplename}.readgroups.bam \
        ID="id" \
        LB="library" \
        PL="illumina" \
        PU="barcode" \
        SM="~{samplename}" \
        CREATE_INDEX=true
      bam=~{samplename}.readgroups.bam
    else
      bam=~{bam}
    fi

    # Generate GVCF depicting single-nucleotide polymorphisms (SNPs) and structural variants (SVs)
    # via local de-novo assembly of haplotypes in variant regions
    gatk \
      --java-options "-Xms~{memory}G -Xmx~{memory}G" \
      HaplotypeCaller \
      -R ${local_ref} \
      -I ${bam} \
      -O ~{samplename}_haplotypecall.g.vcf.gz \
      -ploidy ~{ploidy} \
      --sample-name ~{samplename} \
      ~{"-L " + intervals_file} \
      --tmp-dir . \
      -ERC GVCF

    # CombineGVCF here if multiple read group support is enabled, then feed downstream

    # call GenotypeGVCFs
    # "This tool is able to handle any ploidy (or mix of ploidies) intelligently;
    # there is no need to specify ploidy for non-diploid organisms."
    gatk --java-options "-Xmx~{memory}G" \
      GenotypeGVCFs \
      -R ${local_ref} \
      -V ~{samplename}_haplotypecall.g.vcf.gz \
      ~{"-L " + intervals_file} \
      -O ~{samplename}_genotype.g.vcf.gz
  >>>
  output {
    String gatk_version = read_string("VERSION")
    File gatk_genotype_gvcf = "~{samplename}_genotype.g.vcf.gz"
    File gatk_genotype_gvcf_index = "~{samplename}_genotype.g.vcf.gz.tbi"
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
