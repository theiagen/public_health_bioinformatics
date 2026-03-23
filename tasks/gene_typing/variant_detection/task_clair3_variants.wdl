version 1.0

task clair3_variants {
    input {
    File alignment_bam_file
    File alignment_bam_file_index
    File reference_genome_file
    File reference_genome_file_index
    String sequencing_platform
    String samplename
    String clair3_model = "r1041_e82_400bps_sup_v500"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/clair3-extra-models:1.0.10"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
    Boolean include_all_contigs = true # Haploid calling option, default true
    Boolean enable_haploid_precise = true # Haploid calling option, default true
    Boolean disable_phasing = true # Haploid calling option, default true
    Boolean enable_gvcf = false
    Boolean enable_long_indel = false
    Int variant_quality = 2
    Float snp_min_af = 0.08 # Minimum SNP AF required for a candidate variant.
    Float? indel_min_af # Minimum Indel AF required for a candidate variant. Default dependent on sequencing platform
  }
  String ref_basename = basename(reference_genome_file)
  command <<<
    set -euo pipefail

    # Capture version
    run_clair3.sh --version | head -1 | tee VERSION

    # Set model path
    model_path="/clair3/models/~{clair3_model}"

    # Make sure alignment bam file and fasta index are localized
    echo "Running Clair3 variant calling, aligninment bam file index localized: ~{alignment_bam_file_index}"
    echo "Running Clair3 variant calling, reference genome file index localized: ~{reference_genome_file_index}"

    # Create local fasta & fai copies, Clair3 expects fai & fasta 
    # in working directory, but fai not set explicitly, bam and bai coming from same task
    # so we can assume they are in the same directory
    cp ~{reference_genome_file} ~{ref_basename}
    cp ~{reference_genome_file_index} ~{ref_basename}.fai

    run_clair3.sh \
        --bam_fn=~{alignment_bam_file} \
        --ref_fn=~{ref_basename} \
        --threads=~{cpu} \
        --platform=~{sequencing_platform} \
        --model_path="$model_path" \
        --output=~{samplename} \
        --sample_name=~{samplename} \
        --qual=~{variant_quality} \
        ~{"--snp_min_af=" + snp_min_af} \
        ~{"--indel_min_af=" + indel_min_af} \
        ~{true="--include_all_ctgs" false="" include_all_contigs} \
        ~{true="--haploid_precise" false="" enable_haploid_precise} \
        ~{true="--no_phasing_for_fa" false="" disable_phasing} \
        ~{true="--enable_long_indel" false="" enable_long_indel} \
        ~{true="--gvcf" false="" enable_gvcf}

    gunzip ~{samplename}/merge_output.vcf.gz
    mv "~{samplename}/merge_output.vcf" ~{samplename}_merge_output.vcf

    # If gvcf is enabled, move the gvcf file to the output directory
    if [ "~{enable_gvcf}" == "true" ]; then
      gunzip ~{samplename}/merge_output.gvcf.gz
      mv "~{samplename}/merge_output.gvcf" ~{samplename}_merge_output.gvcf
    fi
  >>>
  output {
    String clair3_version = read_string("VERSION")
    File clair3_variants_vcf = "~{samplename}_merge_output.vcf"
    File? clair3_variants_gvcf = "~{samplename}_merge_output.gvcf"
    String clair3_variants_docker_image = docker
    String clair3_model_used = clair3_model
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 3
  }
}