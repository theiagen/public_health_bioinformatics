version 1.0

task shovill_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "quay.io/staphb/shovill:1.1.0"
    Int min_contig_length = 200
  }
  command <<<
    shovill --version | head -1 | tee VERSION
    shovill \
    --outdir out \
    --R1 ~{read1_cleaned} \
    --R2 ~{read2_cleaned} \
    --minlen ~{min_contig_length}
    mv out/contigs.fa out/~{samplename}_contigs.fasta
    mv out/contigs.gfa out/~{samplename}_contigs.gfa
  >>>
  output {
	  File assembly_fasta = "out/~{samplename}_contigs.fasta"
	  File contigs_gfa = "out/~{samplename}_contigs.gfa"
    String shovill_version = read_string("VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "16 GB"
      cpu: 4
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}

task shovill_se {
  input {
    File read1_cleaned
    String samplename
    String docker = "quay.io/staphb/shovill-se:1.1.0"
    Int min_contig_length = 200
  }
  command <<<
    shovill-se --version | head -1 | tee VERSION
    shovill-se \
    --outdir out \
    --se ~{read1_cleaned} 
    --minlen ~{min_contig_length}
    mv out/contigs.fa out/~{samplename}_contigs.fasta
    mv out/contigs.gfa out/~{samplename}_contigs.gfa
  >>>
  output {
	  File assembly_fasta = "out/~{samplename}_contigs.fasta"
	  File contigs_gfa = "out/~{samplename}_contigs.gfa"
    String shovill_version = read_string("VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "16 GB"
      cpu: 4
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}