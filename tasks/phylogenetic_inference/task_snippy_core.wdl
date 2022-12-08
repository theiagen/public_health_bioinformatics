version 1.0

task snippy_core {
  input {
    Array[File] snippy_outputs
    String treename
    String docker = "staphb/snippy:4.6.0"
    File reference
    File? bed_file
  }
  command <<<
    # version control
   snippy --version | head -1 | tee VERSION
   
   # run snippy core
   snippy-core \
   --prefix ~{treename} \
   ~{'--mask ' + bed_file} \
   --ref ~{reference} \
   ~{sep=" " snippy_outputs}

   # run snippy clean
   snippy-clean_full_aln \
   ~{treename}.full.aln > ~{treename}_snippy_clean_full.aln

  >>>
  output {
   File snippy_core_alignment = "~{treename}.aln"
   File snippy_full_alignment = "~{treename}.full.aln"
   File snippy_full_alignment_clean = "~{treename}.full.aln"
   File snippy_ref = "~{treename}.ref.fa"
   File snippy_core_tab = "~{treename}.tab"
   File snippy_txt = "~{treename}.txt"
   File snippy_vcf = "~{treename}.vcf"
   String snippy_docker_image = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: 4
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
