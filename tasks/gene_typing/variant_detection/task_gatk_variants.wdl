version 1.0

task gatk_variants {
  input {
    String samplename
    File bam
    File bai
    File reference_genome
    Int ploidy = 1 # integer indicating ploidy (N); default to haploid
    File? intervals_file

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/gatk:4.6.2.0-dev"
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
      picard AddOrReplaceReadGroups -Xmx~{memory}G \
        INPUT=~{bam} \
        OUTPUT=~{samplename}.readgroups.bam \
        ID="id" \
        LB="library" \
        PL="illumina" \
        PU="barcode" \
        SM="sample" \
        CREATE_INDEX=true
      bam=~{samplename}.readgroups.bam
      gatk GetSampleName \
        -I ${bam} \
        -O read_groups.txt
    else
      bam=~{bam}
    fi

    # iteratively call HaplotypeCaller on each read_group to generate intermediate GVCF outputs
    mkdir gvcfs/
    count=0
    for read_group in $(cat read_groups.txt | grep -v "\*"); do
      # Generate GVCF depicting single-nucleotide polymorphisms (SNPs) and structural variants (SVs)
      # via local de-novo assembly of haplotypes in variant regions
      gatk \
        --java-options "-Xms~{memory}G -Xmx~{memory}G" \
        HaplotypeCaller \
        -R ${local_ref} \
        -I ${bam} \
        -O gvcfs/~{samplename}_haplotypecall.${count}.g.vcf.gz \
        -ploidy ~{ploidy} \
        --sample-name ${read_group} \
        ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
        --tmp-dir . \
        -ERC GVCF
      count=$((count + 1))
    done

    # CombineGVCFs
    # collect every per-read-group GVCF and pass each as a discrete --variant so all are merged
    vcf_args=""
    for vcf in gvcfs/*.g.vcf.gz; do
      vcf_args="${vcf_args} -V ${vcf}"
    done
    gatk \
      --java-options "-Xmx~{memory}G" \
      CombineGVCFs \
      -R ${local_ref} \
      -O ~{samplename}_haplotypecall.combined.g.vcf.gz \
      ${vcf_args}

    # call GenotypeGVCFs
    # "This tool is able to handle any ploidy (or mix of ploidies) intelligently; 
    # there is no need to specify ploidy for non-diploid organisms."
    gatk --java-options "-Xmx~{memory}G" \
      GenotypeGVCFs \
      -R ${local_ref} \
      -V ~{samplename}_haplotypecall.combined.g.vcf.gz \
      ~{if defined(intervals_file) then "-L ~{intervals_file}" else ""} \
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
