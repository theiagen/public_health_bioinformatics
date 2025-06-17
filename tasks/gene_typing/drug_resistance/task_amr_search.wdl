version 1.0

task amr_search {
  input {
    File input_fasta
    String samplename
    String amr_search_database = "485"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/amrsearch:0.2.1"
    Int cpu = 2
    Int disk_size = 50
    Int memory = 8
  }
  command <<< 
    # Extract base name without path or extension
    input_base=$(basename ~{input_fasta})
    input_base=${input_base%.*}
    echo "DEBUG: input_base = $input_base"

    # Run the tool
    java -jar /paarsnp/paarsnp.jar \
        -i ~{input_fasta} \
        -s ~{amr_search_database}

    # Move the output file from the input directory to the working directory
    mv $(dirname ~{input_fasta})/${input_base}_paarsnp.jsn ./~{samplename}_paarsnp_results.jsn

    # Script housed within the image; https://github.com/theiagen/theiagen_docker_builds/tree/awh-amrsearch-image/amrsearch/0.0.20
    python3 /scripts/parse_amr_json.py \
        ./~{samplename}_paarsnp_results.jsn \
        ~{samplename}
  >>>
  output {
    File json_output = "~{samplename}_paarsnp_results.jsn"
    File output_csv = "~{samplename}_amr_results.csv"
    File output_pdf = "~{samplename}_amr_results.pdf"
    File output_version = "output_amr_version.txt"
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
