version 1.0

task freebayes {
  input {
    String samplename
    File bam
    File bai
    File reference_genome
    Int ploidy = 1 # integer indicating ploidy (N); default to haploid

    File? targets_bed # BED file to limit analysis to specific regions (freebayes -t)

    # native freebayes input filters; left undefined to fall back on freebayes defaults
    Int? min_mapping_quality = 30 # -m: exclude alignments with mapping quality below this (freebayes default: 1)
    Int? min_base_quality = 20 # -q: exclude alleles with supporting base quality below this (freebayes default: 0)
    Float? min_alternate_fraction # -F: min fraction of observations supporting an alternate allele (freebayes default: 0.05)
    Int? min_alternate_count # -C: min count of observations supporting an alternate allele (freebayes default: 2)
    Int? min_coverage # --min-coverage: min coverage required to process a site (freebayes default: 0)

    # gVCF hand-off: freebayes has no native post-call VCF filtering (QUAL/DP/MQ/QD FILTER flags),
    # so emit a bgzipped, tabix-indexed gVCF that can be fed to task_gatk_filter.wdl for that step
    Boolean output_gvcf = false

    String docker = "us-docker.pkg.dev/general-theiagen/staphb/freebayes:1.3.10"
    Int cpu = 8
    Int memory = 32
    Int disk_size = 100
  }
  command <<<
    # fail hard
    set -euo pipefail

    # obtain version (e.g. "version:  v1.3.10" -> "v1.3.10")
    freebayes --version | grep -Po "v[^ ]+$" | tee VERSION

    # localize the reference (symlink in case the source is not writeable) and index it;
    # freebayes needs an accompanying .fai to random-access the reference
    local_ref=$(basename ~{reference_genome})
    ln -s ~{reference_genome} "${local_ref}"
    samtools faidx "${local_ref}"

    # localize the BAM alongside its index so freebayes can find <bam>.bai when a
    # targets BED restricts analysis to specific regions (random access needs the index)
    local_bam=$(basename ~{bam})
    ln -s ~{bam} "${local_bam}"
    ln -s ~{bai} "${local_bam}.bai"

    # shared calling + native input-filter options, assembled once so the VCF and the
    # optional gVCF pass are called identically; undefined optionals render as empty
    # (whitespace-only) tokens and are dropped from the unquoted array
    freebayes_opts=( \
      -f "${local_ref}" \
      ~{if defined(ploidy) then "--ploidy " + ploidy else ""} \
      ~{if defined(min_mapping_quality) then "--min-mapping-quality " + min_mapping_quality else ""} \
      ~{if defined(min_base_quality) then "--min-base-quality " + min_base_quality else ""} \
      ~{if defined(min_alternate_fraction) then "--min-alternate-fraction " + min_alternate_fraction else ""} \
      ~{if defined(min_alternate_count) then "--min-alternate-count " + min_alternate_count else ""} \
      ~{if defined(min_coverage) then "--min-coverage " + min_coverage else ""} \
      ~{if defined(targets_bed) then "--targets " + targets_bed else ""} \
    )

    # primary variant call: standard VCF output, then bgzip + tabix index
    freebayes "${freebayes_opts[@]}" "${local_bam}" > ~{samplename}.freebayes.vcf
    bgzip ~{samplename}.freebayes.vcf
    tabix -p vcf ~{samplename}.freebayes.vcf.gz

    # optional gVCF pass (identical calling options as the VCF pass), emitting a
    # bgzipped, tabix-indexed .g.vcf.gz + .tbi pair suitable for GATK (task_gatk_filter).
    # freebayes --gvcf emits each contig's trailing reference block as a malformed record
    # with a missing reference allele (REF ".") carrying the previous contig's boundary
    # coordinates. GATK/htsjdk rejects "REF ." ("reference allele cannot be missing") and
    # those same records also sit out of positional order (breaking tabix). They are
    # non-variant blocks, so drop them (awk), then restore positional order.
    if [ "~{output_gvcf}" == "true" ]; then
      # vcfstreamsort -a loads the whole gvcf in memory, which ought to be fine for non-plant/animal genomes
      freebayes "${freebayes_opts[@]}" --gvcf "${local_bam}" \
        | awk -F'\t' '/^#/ || $4 != "."' \
        | vcfstreamsort -a > ~{samplename}.freebayes.g.vcf
      bgzip ~{samplename}.freebayes.g.vcf
      tabix -p vcf ~{samplename}.freebayes.g.vcf.gz
    fi
  >>>
  output {
    String freebayes_version = read_string("VERSION")
    String freebayes_docker_image = docker
    File freebayes_vcf = "~{samplename}.freebayes.vcf.gz"
    File freebayes_vcf_index = "~{samplename}.freebayes.vcf.gz.tbi"
    File? freebayes_gvcf = "~{samplename}.freebayes.g.vcf.gz"
    File? freebayes_gvcf_index = "~{samplename}.freebayes.g.vcf.gz.tbi"
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
