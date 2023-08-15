version 1.0

task czgenepi_wrangling {
  input {
    Array[String] sample_names
    Array[File] assembly_fasta
    Array[String] collection_date
    Array[String] private_id

    # collection location
    Array[String] continent
    Array[String] country
    Array[String] state
    Array[String]? county

    # optional inputs
    Array[String]? gisaid_virus_name
    Array[String]? sequencing_date
    Array[String]? sample_is_private

    # runtime
    Int disk_size = 100
  }
  command <<<
  >>>
  output {
    File concatenated_fasta = "file"
    File concatenated_metadata = "file"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}