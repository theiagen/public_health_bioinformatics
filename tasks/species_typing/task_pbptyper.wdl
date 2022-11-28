version 1.0

task pbptyper {
  meta {
    description: " In silico Penicillin Binding Protein (PBP) typer for Streptococcus pneumoniae assemblies. https://github.com/rpetit3/pbptyper"
  }
  input {
    File assembly # An assembly in FASTA format (compressed with gzip, or uncompressed) to predict the PBP type on.
    String samplename
    String? db # A path to a directory containing FASTA files for 1A, 2B, and 2X proteins. In most cases using the default value will be all that is needed.
    Int min_pident = 95 # Minimum percent identity to count a hit [default: 95]
    Int min_coverage = 95 # Minimum percent coverage to count a hit [default: 95]  
    String docker = "staphb/pbptyper:1.0.4"
    Int cpus = 4

  }
  command <<<
    # get version information
    pbptyper --version | sed 's/pbptyper, //' | tee VERSION
    
    # run pbptyper
    pbptyper \
      --assembly ~{assembly} \
      ~{'--db ' + db} \
      ~{'--min_pident ' + min_pident} \
      ~{'--min_coverage ' + min_coverage} \
      --prefix "~{samplename}" \
      --outdir ./ 

    # parse output tsv for pbptype
    cut -f 2 ~{samplename}.tsv | tail -n 1 > pbptype.txt

  >>>
  output {
    String pbptyper_predicted_1A_2B_2X = read_string("pbptype.txt")
    File pbptyper_pbptype_predicted_tsv = "~{samplename}.tsv" # A tab-delimited file with the predicted PBP type
    File pbptyper_pbptype_1A_tsv = "~{samplename}-1A.tblastn.tsv" # A tab-delimited file of all blast hits against 1A
    File pbptyper_pbptype_2B_tsv = "~{samplename}-2B.tblastn.tsv" # A tab-delimited file of all blast hits against 2B
    File pbptyper_pbptype_2X_tsv = "~{samplename}-2X.tblastn.tsv" # A tab-delimited file of all blast hits against 2X
    String pbptyper_version = read_string("VERSION")
    String pbptyper_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: "16 GB"
    cpu: cpus
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
