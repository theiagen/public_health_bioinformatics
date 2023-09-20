version 1.0

task prep_theiacov_fasta_metadata {
  input {
    File assembly
    String seq_method
    String assembly_method
    String organism
    String flu_segment

    Int disk_size = 10
  }
  command <<<
    # Set strain name by assembly header
    assembly_header=$(grep -e ">" ~{assembly} | sed 's/\s.*$//' |  sed 's/>//g' )

    sequencing_method=""
    assembly_method=""
    influenza_segment="influenza_segment"

    # if seq_method defined, add to metadata
    if [[ "~{seq_method}" ]]; then
      sequencing_method="seq_method"
    fi
    # if assembly_method defined, add to metadata
    if [[ "~{assembly_method}" ]]; then
      assembly_method="assembly_method"
    fi
    # if organism defined, add to metadata
    if [[ "~{organism}" == "sars-cov-2" ]]; then
      organism="sars-cov-2"
    else
      organism="~{organism}"
    fi
    # if flu segment defined, add to metadata
    if [[ "~{flu_segment}" == "HA" ]]; then
      flu_segment="HA"
    else
      flu_segment="~{flu_segment}"
    fi

    # write everything to a file
    echo -e "assembly\t${sequencing_method}\t${assembly_method}\torganism\t${influenza_segment}" > theaicov_fasta_set_metadata.tsv
    echo -e "\"${assembly_header}\"\t\"~{seq_method}\"\t\"~{assembly_method}\"\t\"${organism}\"\t\"${flu_segment}\"" >> theaicov_fasta_set_metadata.tsv
  >>>
  output {
    File theaicov_fasta_set_metadata = "theaicov_fasta_set_metadata.tsv"
  }
  runtime {
      docker: "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
      memory: "3 GB"
      cpu: 1
      disks: "local-disk ~{disk_size} SSD"
      disk: disk_size + " GB"
      preemptible: 0
      maxRetries: 3
  }
}