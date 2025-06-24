version 1.0

task bcftools_consensus {
  input{
    File reference_fasta
    File input_vcf
    String samplename
    Int min_depth = 0
    Float min_freq = 0.0
    Int disk_size = 100
    Int cpu = 2
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/bcftools:1.20"
  }
  command <<<
    set -euo pipefail

    # get version
    bcftools --version | head -1 | tee VERSION

    # remove low coverage variants (this is where you would remove quality too, e.g. & MIN(FMT/GQ)>MIN_QUAL)
    echo "DEBUG: Filtering variants with low coverage: ~{min_depth} and low frequency: ~{min_freq}"
    bcftools view -i 'MIN(FMT/DP)>=~{min_depth} && MIN(FMT/AF)>=~{min_freq}' \
      ~{input_vcf} \
      > ~{samplename}_cov_filtered_prefinal.vcf

    # reproduce artic behavior for left-aligning and normalizing indels
    echo "Left-normalizing indels"
    bcftools norm \
      ~{samplename}_cov_filtered_prefinal.vcf \
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
      --output ~{samplename}_temp.fasta

    # prepend samplename to fasta headers
    sed -e 's/^>/>'~{samplename}'_/' ~{samplename}_temp.fasta > ~{samplename}_consensus.fasta

    # decompress the vcf file
    gunzip ~{samplename}_filtered.vcf.gz
  >>>
  output {
    File assembly_fasta = "~{samplename}_consensus.fasta"
    File bcftools_filtered_vcf = "~{samplename}_filtered.vcf"
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