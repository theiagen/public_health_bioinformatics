version 1.0

task trim_genbank_fastas {
  input {
    File genbank_untrimmed_fasta
    String output_name
    Int min_length = 50
    Int max_length = 30000
    Int disk_size = 100
    Int memory = 2
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/vadr:1.3"
  }
  command <<<
    # remove terminal ambiguous nucleotides
    /opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
      ~{genbank_untrimmed_fasta} \
      --minlen ~{min_length} \
      --maxlen ~{max_length} \
      > ~{output_name}_genbank.fasta
  >>>
  output {
    File genbank_fasta = "~{output_name}_genbank.fasta"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}


## I think this works, but honestly not sure.
task table2asn {
  input {
    File authors_sbt # have users provide the .sbt file for MPXV submission-- it can be created here: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
    File bankit_fasta
    File bankit_metadata
    String output_name
    Int disk_size = 100
    Int memory = 1
    Int cpu = 1
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/ncbi-table2asn:1.26.678"
  }
  command <<<
    # using this echo statement so the fasta file doesn't have a wiggly line
    echo "~{bankit_fasta} file needs to be localized for the program to access"

    # rename authors_sbt to contain output_name so table2asn can find it
    # had issues with device busy, making softlinks for all
    ln -s ~{authors_sbt} ~{output_name}.sbt
    ln -s ~{bankit_fasta} ~{output_name}.fsa
    ln -s ~{bankit_metadata} ~{output_name}.src

    # convert the data into a sqn file so it can be emailed to NCBI
    table2asn -t ~{output_name}.sbt \
      -src-file ~{output_name}.src \
      -indir . \
      -a s # inputting a set of fasta data
    
  >>>
  output {
    File sqn_file = "~{output_name}.sqn"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
    continueOnReturnCode: [0, 2]
  }
}