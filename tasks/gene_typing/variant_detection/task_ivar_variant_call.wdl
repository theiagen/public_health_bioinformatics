version 1.0

task variant_call {
  input {
    File bamfile
    String samplename
    File? reference_genome
    File? reference_gff 
    Boolean count_orphans = true
    Int max_depth = "600000"
    Boolean disable_baq = true
    Int min_bq = "0"
    Int min_qual = "20"
    Float? variant_min_freq 
    Int? variant_min_depth 
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ivar:1.3.1-titan"
  }
  command <<<
    set -euo pipefail
    # version control
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
    
    # set reference gff
    if [[ ! -z "~{reference_gff}" ]]; then
      echo "User reference identified; ~{reference_genome} will be utilized for alignement"
      ref_gff="~{reference_gff}"
      # move to primer_schemes dir; bwa fails if reference file not in this location
    else
      ref_gff="/reference/GCF_009858895.2_ASM985889v3_genomic.gff"  
    fi
    
    # call variants
    samtools mpileup \
    ~{true = "-A" false = "" count_orphans} \
    -d ~{max_depth} \
    ~{true = "-B" false = "" disable_baq} \
    -Q ~{min_bq} \
    --reference ${ref_genome} \
    ~{bamfile} | \
    ivar variants \
    -p ~{samplename}.variants \
    -q ~{min_qual} \
    -t ~{variant_min_freq} \
    -m ~{variant_min_depth} \
    -r ${ref_genome} \
    -g ${ref_gff}

    # Convert TSV to VCF
    ivar_variants_to_vcf.py ~{samplename}.variants.tsv ~{samplename}.variants.vcf

    # Variant calculations

    # keep only unique nucleotide substitutions:
    # the same nucleotide substitution can be listed in multiple rows if it is within 
    # more than one coding region, so remove columns with coding region information 
    # then keep only unique rows
    cut -f 1-14 ~{samplename}.variants.tsv | awk '!seen[$0]++ {print}' > unique_variants.tsv

    # filter variants that pass fisher's exact test
    grep "TRUE" unique_variants.tsv > passed_variants.tsv || touch passed_variants.tsv

    # calculate total number of variants
    variants_num=$(wc -l < passed_variants.tsv || echo 0)
    echo $variants_num | tee VARIANT_NUM

    # calculate proportion of variants with allele frequencies between 0.6 and 0.9
    # find number of variants at intermediate frequencies 
    awk -F "\t" '{ if(($11 >= 0.6) && ($11 <= 0.9)) {print }}' passed_variants.tsv > intermediate_variants.tsv
    intermediates_num=$(wc -l < intermediate_variants.tsv || echo 0)

    # if number of total variants is not zero, divide number of intermediate variants by total number of variants
    if [[ "$variants_num" -eq "0" ]]; then
      echo "Not computed: no variants detected" > PROPORTION_INTERMEDIATE
    else
      echo $intermediates_num $variants_num | awk '{ print $1/$2 }' > PROPORTION_INTERMEDIATE
    fi
  >>>
  output {
    Int variant_num = read_string("VARIANT_NUM")
    String variant_proportion_intermediate = read_string("PROPORTION_INTERMEDIATE")
    File sample_variants_tsv = "~{samplename}.variants.tsv"
    File sample_variants_vcf = "~{samplename}.variants.vcf"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
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
