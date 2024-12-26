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
    Boolean include_all_contigs = true # Haploid calling option, default true
    Boolean enable_haploid_precise = true # Haploid calling option, default true
    Boolean disable_phasing = true # Haploid calling option, default true
    Boolean enable_gvcf = false
    Boolean enable_long_indel = false
    Int variant_quality = 2
  }
  String bam_basename = basename(alignment_bam_file)
  String ref_basename = basename(reference_genome_file)
  command <<<
    set -euo pipefail

    model_path="/clair3/models/~{clair3_model}"

    echo "Running Clair3 variant calling, aligninment bam file index localized: ~{alignment_bam_file_index}"
    echo "Running Clair3 variant calling, reference genome file index localized: ~{reference_genome_file_index}"

    # Create local copies with correct names, Clair3 expects bai and fai in working directory, but not set explicitly
    cp ~{alignment_bam_file} ~{bam_basename}
    cp ~{alignment_bam_file_index} ~{bam_basename}.bai
    cp ~{reference_genome_file} ~{ref_basename}
    cp ~{reference_genome_file_index} ~{ref_basename}.fai

    run_clair3.sh \
        --bam_fn=~{bam_basename} \
        --ref_fn=~{ref_basename} \
        --threads=~{cpu} \
        --platform=~{sequencing_platform} \
        --model_path="$model_path" \
        --output=~{samplename} \
        --sample_name=~{samplename} \
        --qual=~{variant_quality} \
        ~{true="--include_all_ctgs" false="" include_all_contigs} \
        ~{true="--haploid_precise" false="" enable_haploid_precise} \
        ~{true="--no_phasing_for_fa" false="" disable_phasing} \
        ~{true="--enable_long_indel" false="" enable_long_indel} \
        ~{true="--gvcf" false="" enable_gvcf}

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