version 1.0

task parse_paf {
  meta {
    description: "Parse minimap2 PAF file and return aligned contigs in FASTA format"
  }
  input {
    File paf
    File assembly
    String samplename
    String docker = "staphb/seqkit:2.3.1"
    Int disk_size = 100
    Int cpu = 2
    Int mem = 8
  }
  command <<<
    # retrieve contig name
    cut -f1 "~{paf}" > contig_names.txt
    
    # extract mapped contigs in FASTA format
    seqkit grep -f contig_names.txt "~{assembly}" > "~{samplename}"_mapped_contigs.fasta

  >>>
  output {
    File parse_paf_contigs = "~{samplename}_mapped_contigs .fasta"
  }
  runtime {
    docker: "~{docker}"
    memory: mem + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}