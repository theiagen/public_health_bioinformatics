version 1.0

task dragonflye {
  input {
    File reads
  }
  command <<<
    # do dragonflye woo
  >>>
  output {
    File assembly_fasta 
  }
  runtime {

  }
}