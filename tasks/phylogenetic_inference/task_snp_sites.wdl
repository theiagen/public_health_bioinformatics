version 1.0

task snp_sites {
  input {
    File alignment
    String samplename
    String docker = "staphb/snp-sites:2.5.1"
    Int cpu = 4
    Int disk_size = 100
  }
  command <<<
    # date and version control
    date | tee DATE
    snp-sites -V | tee VERSION

    snp-sites \
    -cvpm \
    -o ~{samplename} \
    ~{alignment}
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    String snpsites_docker_image = docker
    File snpsites_fasta = "~{samplename}.snp_sites.aln"
    File snpsites_vcf = "~{samplename}.vcf"
    File snpsites_phylip = "~{samplename}.phylip"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 1
  }
}
