version 1.0

task kmerfinder_bacteria {
  input {
    File assembly
    String samplename
    File kmerfinder_db = "gs://theiagen-public-resources-rp/reference_data/databases/kmerfinder/kmerfinder_bacteria_20230911.tar.gz"
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/kmerfinder:3.0.2--hdfd78af_0"
    Int memory = 32
    Int cpu = 4
    Int disk_size = 100
    String kmerfinder_args = ""
  }
  command <<<
    # Decompress the kmerfinder bacterial database
    mkdir db
    tar -C ./db/ -xzvf ~{kmerfinder_db}  

    # Run kmerfinder
    kmerfinder.py \
        -db ./db/bacteria/bacteria.ATG \
        -tax ./db/bacteria/bacteria.tax \
        -i ~{assembly} \
        -o ~{samplename} \
        ~{kmerfinder_args} 

    # parse outputs
    if [ ! -f ~{samplename}/results.txt ]; then
      PF="No hit detected in database"
      QC="No hit detected in database"
      TC="No hit detected in database"
    else
      PF="$(cat ~{samplename}/results.txt | head -n 2 | tail -n 1 | cut -f 19)"
      QC="$(cat ~{samplename}/results.txt | head -n 2 | tail -n 1 | cut -f 6)"
      TC="$(cat ~{samplename}/results.txt | head -n 2 | tail -n 1 | cut -f 7)"
        # String is empty or just contains the header
        if [ "$PF" == "" ] || [ "$PF" == "Species" ]; then
          PF="No hit detected in database"
          QC="No hit detected in database"
          TC="No hit detected in database"
        fi
      mv -v ~{samplename}/results.txt ~{samplename}_kmerfinder.tsv
    fi
    echo $PF | tee TOP_HIT
    echo $QC | tee QC_METRIC
    echo $TC | tee TEMPLATE_COVERAGE

    # extract database name
    DB=$(basename ~{kmerfinder_db} | sed 's/\.tar\.gz$//')
    echo $DB | tee DATABASE
  >>>
  output {
    String kmerfinder_docker = docker
    File? kmerfinder_results_tsv = "~{samplename}_kmerfinder.tsv"
    String kmerfinder_top_hit = read_string("TOP_HIT")
    String kmerfinder_query_coverage = read_string("QC_METRIC")
    String kmerfinder_template_coverage = read_string("TEMPLATE_COVERAGE")
    String kmerfinder_database = read_string("DATABASE")
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk ~{disk_size} SSD"
    preemptible: 0
  }
}