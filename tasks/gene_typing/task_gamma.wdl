version 1.0

task gamma_one_sample {
  input {
    File assembly_fasta
    String samplename
    String docker = "quay.io/biocontainers/gamma:1.4--hdfd78af_0"
    File gamma_database
  }
  String database_name = basename(gamma_database)
  command <<<
    GAMMA.py ~{assembly_fasta} ~{gamma_database} ~{samplename}
    
    mv ~{samplename}.gamma ~{samplename}_gamma_report.tsv
    
  >>>
  output {
    File gamma_results = "~{samplename}_gamma_report.tsv"
    String gamma_database_version = database_name
    String gamma_docker = docker
  }
  runtime {
    memory: "8 GB"
    cpu: 4
    docker: "~{docker}"
    disks: "local-disk 100 HDD"
  }
}