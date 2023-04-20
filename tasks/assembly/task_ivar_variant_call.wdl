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
  }
  command <<<
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
    -t ~{min_freq} \
    -m ~{variant_min_depth} \
    -r ${ref_genome} \
    -g ${ref_gff}

    # Convert TSV to VCF
    ivar_variants_to_vcf.py ~{samplename}.variants.tsv ~{samplename}.variants.vcf

    variants_num=$(grep "TRUE" ~{samplename}.variants.tsv | wc -l)
    if [ -z "$variants_num" ] ; then variants_num="0" ; fi
    echo $variants_num | tee VARIANT_NUM
  >>>
  output {
    Int variant_num = read_string("VARIANT_NUM")
    File sample_variants_tsv = "~{samplename}.variants.tsv"
    File sample_variants_vcf = "~{samplename}.variants.vcf"
    String ivar_version = read_string("IVAR_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
    String pipeline_date = read_string("DATE")
  }
  runtime {
    docker: "quay.io/staphb/ivar:1.3.1-titan"
    memory: "8 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}