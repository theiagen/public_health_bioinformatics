version 1.0

task medaka_polish {
    input {
      File assembly_fasta
      String samplename
      File read1
      String medaka_model
      Int cpu = 4
      Int memory = 16
      Int disk_size = 100
      String docker = "staphb/medaka:2.0.1"
    }
    command <<< 
        medaka --version | tee VERSION
        medaka_consensus \
          -i ~{read1} \
          -d ~{assembly_fasta} \
          -m ~{medaka_model} \
          -o . \
          --threads ~{cpu}
        # Rename output files with sample name
        mv consensus.fasta ~{samplename}.polished.fasta
        mv consensus.vcf.gz ~{samplename}.polished.vcf.gz
    >>>
    output {
      File polished_fasta = "~{samplename}.polished.fasta"
      File polished_vcf = "~{samplename}.polished.vcf.gz"
      String medaka_version = read_string("VERSION")
    }
    runtime {
      docker: "~{docker}"
      cpu: cpu
      memory: "~{memory} GB"
      disks: "local-disk " + disk_size + " HDD"
      disk: disk_size + " GB"
      maxRetries: 3
      preemptible: 0
    }
}
