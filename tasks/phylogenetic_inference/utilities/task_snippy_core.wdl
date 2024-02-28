version 1.0

task snippy_core {
  input {
    Array[File] snippy_variants_outdir_tarball
    Array[String] samplenames
    String tree_name
    File reference_genome_file
    File? bed_file

    String docker = "us-docker.pkg.dev/general-theiagen/staphb/snippy:4.6.0"
    Int disk_size = 100
    Int cpu = 8
    Int memory = 16
  }
  command <<<
    # version control
    snippy --version | head -1 | tee VERSION

    tarball_array=(~{sep=" " snippy_variants_outdir_tarball})
    samplename_array=(~{sep=" " samplenames})

    # iteratively untar
    for i in ${tarball_array[@]}; do tar -xf $i; done

    # run snippy core
    snippy-core \
    --prefix ~{tree_name} \
    ~{'--mask ' + bed_file} \
    --ref ~{reference_genome_file} \
    "${samplename_array[@]}"

    # run snippy clean
    snippy-clean_full_aln \
      ~{tree_name}.full.aln > ~{tree_name}_snippy_clean_full.aln

    # commenting out since the _core.aln should not be confused with the core genome alignment created later in snippy_tree workflow
    #mv ~{tree_name}.aln ~{tree_name}_core.aln
    mv ~{tree_name}.full.aln ~{tree_name}_full.aln
    mv ~{tree_name}.tab ~{tree_name}_all_snps.tsv
    mv ~{tree_name}.txt ~{tree_name}_snps_summary.txt
  >>>
  output {
    String snippy_core_version = read_string("VERSION")
    # commenting out for same reason as above
    #File snippy_core_alignment = "~{tree_name}_core.aln"
    File snippy_full_alignment = "~{tree_name}_full.aln"
    File snippy_full_alignment_clean = "~{tree_name}_snippy_clean_full.aln"
    File snippy_ref = "~{tree_name}.ref.fa"
    File snippy_core_tab = "~{tree_name}_all_snps.tsv"
    File snippy_txt = "~{tree_name}_snps_summary.txt"
    File snippy_vcf = "~{tree_name}.vcf"
    String snippy_core_docker_image = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}