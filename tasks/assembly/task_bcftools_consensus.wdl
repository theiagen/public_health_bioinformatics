version 1.0

task bcftools_consensus {
  input{
    File reference_fasta
    File input_vcf
    String samplename
    Int min_depth = 0
    Float min_freq = 0.0
    Int disk_size = 100
    Boolean decompress = false
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bcftools:1.20"
  }
  command <<<
    set -euo pipefail

    # get version
    bcftools --version | head -1 | tee VERSION

    # make sure input vcf is gzipped
    if [[ ! "~{input_vcf}" == *.gz ]]; then
      ext=".gz"
      gzip ~{input_vcf}
    else
      ext=""
    fi

    # remove low coverage variants (this is where you would remove quality too, e.g. & MIN(FMT/GQ)>MIN_QUAL)
    echo "DEBUG: Filtering variants with low coverage: ~{min_depth} and low frequency: ~{min_freq}"
    bcftools view -i 'MIN(FMT/DP)>~{min_depth} && MIN(FMT/AF)>~{min_freq}' \
      ~{input_vcf}${ext} \
      > ~{samplename}_cov_filtered_prefinal.vcf.gz

    # reproduce artic behavior for left-aligning and normalizing indels
    echo "Left-normalizing indels"
    bcftools norm \
      ~{samplename}_cov_filtered_prefinal.vcf.gz \
      --check-ref wx \
      --fasta-ref ~{reference_fasta} \
      --output-type z \
      --output ~{samplename}_filtered.vcf.gz

    # index the vcf file
    echo "Indexing VCF"
    bcftools index ~{samplename}_filtered.vcf.gz

    # create the consensus fasta file
    echo "Generating consensus sequence"
    bcftools consensus \
      ~{samplename}_filtered.vcf.gz \
      --fasta-ref ~{reference_fasta} \
      --output ~{samplename}_consensus.fasta

    if [ "~{decompress}" == "true" ]; then
      echo "DEBUG: Decompressing VCF"
      gunzip ~{samplename}_filtered.vcf.gz
    fi
  >>>
  output {
    File assembly_fasta = "~{samplename}_consensus.fasta"
    File bcftools_filtered_vcf = "~{samplename}_filtered.vcf~{true = "" false = ".gz" decompress}"
    String bcftools_version = read_string("VERSION")
    String bcftools_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
    preemptible: 0
  }
}