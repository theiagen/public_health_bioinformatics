version 1.0

task ksnp3 {
  input {
    File alignment
    File reference
    String sample_name
    String model = "ont"
    String docker = "us-docker.pkg.dev/general-theiagen/hkubal/clair3:v1.0.6"
    Int memory = 8
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    model_path="/opt/models/${model_name}"
    # hifi
    # hifi_revio
    # hifi_sequel2 -> hifi
    # ilmn
    # ont -> r941_prom_hac_g360+g422
    # ont_guppy2 -> r941_prom_hac_g238
    # ont_guppy5 -> r941_prom_sup_g5014
    # r941_prom_hac_g238
    # r941_prom_hac_g360+g422
    # r941_prom_sup_g5014
    tmpoutdir=$(mktemp -d)

    run_clair3.sh \
        --bam_fn=~{alignment} \
        --ref_fn=~{reference} \
        --threads=~{cpu} \
        --platform="ont" \
        --model_path="$model_path" \
        --output="$tmpoutdir" \
        --sample_name=~{sample_name} \
        --include_all_ctgs \
        --haploid_precise \
        --no_phasing_for_fa \
        --enable_long_indel

    mv "${tmpoutdir}/merge_output.vcf.gz" ~{sample_name}.vcf.gz
  >>>
  output {
    File clair3_vcf = "~{sample_name}.vcf.gz"
    String ksnp3_docker_image = docker
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
