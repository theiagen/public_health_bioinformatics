version 1.0

task bcftools_consensus {
  input{
    File reference_fasta
    File input_vcf
    String samplename
    Int min_depth = 0
    Int disk_size = 100
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

    # remove low coverage variants (this is where you would remove quality too, e.g. & MIN(FMT/GQ > MIN_QUAL))
    bcftools view -i 'MIN(FMT/DP>~{min_depth})' \
      ~{input_vcf}${ext} \
      > ~{samplename}_cov_filtered.vcf.gz

    # reproduce artic behavior for left-aligning and normalizing indels
    bcftools norm \
      ~{samplename}_cov_filtered.vcf.gz \
      --check-ref wx \
      --fasta-ref ~{reference_fasta} \
      --output-type z \
      --output ~{samplename}_norm.vcf.gz

    # index the vcf file
    bcftools index ~{samplename}_norm.vcf.gz

    # create the consensus fasta file
    bcftools consensus \
      ~{samplename}_norm.vcf.gz \
      --fasta-ref ~{reference_fasta} \
      --output ~{samplename}_consensus.fasta
  >>>
  output {
    File assembly_fasta = "~{samplename}_consensus.fasta"
    File bcftools_norm_vcf = "~{samplename}_norm.vcf.gz"
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