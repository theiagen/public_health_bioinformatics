version 1.0

task amr_search {
  input {
    File input_fasta
    String samplename
    String database = "485"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.1.0"
    Int cpu = 2
    Int disk_size = 50
    Int memory = 8
  }
  command <<< 
    # Extract base name without path or extension
    input_base=$(basename ~{input_fasta} .fasta)

    # Run the tool
    java -jar /paarsnp/paarsnp.jar \
        -i ~{input_fasta} \
        -s ~{database}

    # Move the output file from the input directory to the working directory
    mv $(dirname ~{input_fasta})/${input_base}_paarsnp.jsn ./~{samplename}_paarsnp_results.jsn
  >>>
  output {
    File json_output = "~{samplename}_paarsnp_results.jsn"
    String amr_search_docker = docker
  }

  runtime {
    memory: "~{memory} GB"
    cpu: cpu
    docker: docker
    disks: "local-disk " + disk_size + " SSD"
    maxRetries: 1
  }
}
