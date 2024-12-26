version 1.0

task clair3_variants {
    input {
    File alignment_bam_file
    File alignment_bam_file_index
    File reference_genome_file
    File reference_genome_file_index
    String sequencing_platform
    String samplename
    String clair3_model = "r941_prom_hac_g360+g422"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/clair3:1.0.10"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
    Boolean enable_gvcf = false
    Boolean enable_long_indel = false
    Int variant_quality = 2
  }
  command <<<
    set -euo pipefail

    model_path="/clair3/models/~{clair3_model}"

    echo "Running Clair3 variant calling, aligninment bam file index localized: ~{alignment_bam_file_index}"
    echo "Running Clair3 variant calling, reference genome file index localized: ~{reference_genome_file_index}"

    run_clair3.sh \
        --bam_fn=~{alignment_bam_file} \
        --ref_fn=~{reference_genome_file} \
        --threads=~{cpu} \
        --platform=~{sequencing_platform} \
        --model_path="$model_path" \
        --output=~{samplename} \
        --sample_name=~{samplename} \
        --qual=~{variant_quality} \
        --include_all_ctgs \
        --haploid_precise \
        --no_phasing_for_fa \
        ~{if enable_long_indel then "--enable_long_indel" else ""} \
        ~{if enable_gvcf then "--gvcf" else ""}

    mv "~{samplename}/merge_output.vcf.gz" ~{samplename}_merge_output.vcf.gz
    mv "~{samplename}/pileup.vcf.gz" ~{samplename}_pileup.vcf.gz
    mv "~{samplename}/full_alignment.vcf.gz" ~{samplename}_full_alignment.vcf.gz

    # If gvcf is enabled, move the gvcf file to the output directory
    if [ "~{enable_gvcf}" == "true" ]; then
        mv "~{samplename}/merge_output.gvcf.gz" ~{samplename}_merge_output.gvcf.gz
    fi
  >>>
  output {
    File clair3_variants_final_vcf = "~{samplename}_merge_output.vcf.gz"
    File clair3_variants_full_alignment_vcf = "~{samplename}_full_alignment.vcf.gz"
    File clair3_variants_pileup_vcf = "~{samplename}_pileup.vcf.gz"
    File? clair3_variants_gvcf = "~{samplename}_merge_output.gvcf.gz"
    String clair3_variants_docker_image = docker
    String clair3_model_used = clair3_model
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
    maxRetries: 0
  }
}