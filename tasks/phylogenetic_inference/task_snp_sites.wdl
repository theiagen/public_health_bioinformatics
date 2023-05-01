version 1.0

task snp_sites {
  input {
    File msa_fasta
    String output_name
    Boolean allow_wildcard_bases = true
    Boolean output_vcf = true
    Boolean output_phylip = false
    Boolean output_multifasta = false
    Boolean output_pseudo_ref = false
    Boolean output_monomorphic = false
    String docker = "quay.io/staphb/snp-sites:2.5.1"
    Int disk_size = 100
    Int cpus = 1
    Int memory = 4
  }
  command <<< 
    snp-sites -V > VERSION

    # Usage: snp-sites [-mvph] [-o output_filename] <file>
    #  -r     output internal pseudo reference sequence
    #  -m     output a multi fasta alignment file
    #  -v     output a VCF file (default)
    #  -p     output a phylip file
    #  -c     only output columns containing exclusively ACGT (false by default)
    #  -b     output monomorphic sites, used for BEAST
    #  <file> input alignment file which can optionally be gzipped

    snp-sites \
      ~{true="-v" false="" output_vcf} \
      ~{true="-p" false="" output_phylip} \
      ~{true="-m" false="" output_multifasta} \
      ~{true="-r" false="" output_pseudo_ref} \
      ~{true="" false="-c" allow_wildcard_bases} \
      ~{true="-b" false="" output_monomorphic} \
      -o "~{output_name}" "~{msa_fasta}"
  >>>
  output {
    File? snp_sites_vcf = "~{output_name}.vcf"
    File? snp_sites_phylip = "~{output_name}.phylip"
    File? snp_sites_multifasta = "~{output_name}.snp_sites.aln"
    String snp_sites_version = read_string("VERSION")
    String snp_sites_docker = docker
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu :  cpus
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x4"
    maxRetries: 3
  }
}