version 1.0

workflow de_novo_assembly {

  input {
    String    SRR
    File      read1
    File      read2
  }

  call seqyclean {
    input:
      samplename=SRR,
      read1=read1,
      read2=read2
  }

  call shovill {
    input:
      samplename=SRR,
      read1_cleaned=seqyclean.read1_cleaned,
      read2_cleaned=seqyclean.read2_cleaned
  }

  output {
    File    read1_cleaned     =seqyclean.read1_cleaned
    File    read2_cleaned     =seqyclean.read2_cleaned
    File    contigs_fasta     =shovill.contigs_fasta
    File    contigs_gfa       =shovill.contigs_gfa
  }
}

task seqyclean {

  input {
    File        read1
    File        read2
    String      samplename
    File?       adapters
    Int?        seqyclean_minlen=25
    String?     seqyclean_qual="20 20"
    Boolean?    compress=true
    Boolean?    seqyclean_dup=false
    Boolean?    seqyclean_no_adapter_trim=false
  }

  command {
    seqyclean --version | head -1 | tee VERSION
    seqyclean \
    ${'-minlen ' + seqyclean_minlen} \
    ${'-qual ' + seqyclean_qual} \
    ${'-c ' + adapters} \
    ${true="-dup" false="" seqyclean_dup} \
    ${true="-no_adapter_trim " false="" seqyclean_no_adapter_trim} \
    ${true="-gz " false="" compress} \
    ${'-1 ' + read1} \
    ${'-2 ' + read2} \
    ${'-o ' + samplename}
  }

  output {
	  File       read1_cleaned = "${samplename}_PE1.fastq.gz"
	  File       read2_cleaned = "${samplename}_PE2.fastq.gz"
    String     seqyclean_version = read_string("VERSION")
  }

  runtime {
      docker:       "quay.io/staphb/seqyclean:1.10.09"
      memory:       "8 GB"
      cpu:          2
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}

task shovill {

  input {
    File        read1_cleaned
    File        read2_cleaned
    String      samplename
  }

  command {
    shovill --version | head -1 | tee VERSION
    shovill \
    --outdir out \
    --R1 ${read1_cleaned} \
    --R2 ${read2_cleaned}
    mv out/contigs.fa out/${samplename}_contigs.fasta
    mv out/contigs.gfa out/${samplename}_contigs.gfa
  }

  output {
	  File       contigs_fasta = "out/${samplename}_contigs.fasta"
	  File       contigs_gfa = "out/${samplename}_contigs.gfa"
    String     shovill_version = read_string("VERSION")
  }

  runtime {
      docker:       "quay.io/staphb/shovill:1.1.0"
      memory:       "16 GB"
      cpu:          4
      disks:        "local-disk 100 SSD"
      preemptible:  0
  }
}
