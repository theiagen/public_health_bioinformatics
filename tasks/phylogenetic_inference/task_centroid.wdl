version 1.0

task centroid {
  input {
    Array[File] assembly_fastas
  }
  command <<<
    python3 centroid.py .
  >>>
  output {
    String centroid_genome = read_string("centroid_out.txt")
  }
  runtime {
    docker: "curtis-docker"
    
  }
}