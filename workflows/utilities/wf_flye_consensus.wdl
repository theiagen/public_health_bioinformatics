version 1.0

workflow flye_consensus {
  meta {
    description: "This workflow runs flye to generate a consensus genome assembly from long reads."
  }
  input {
    File read1
    Int genome_size
    Int threads
  }
}