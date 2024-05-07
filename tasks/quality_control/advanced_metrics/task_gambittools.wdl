version 1.0

task gambit_core_check {
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/internal/gambittools:1.0.0"
    File gambit_db_genomes = "gs://gambit-databases-rp/2.0.0/gambit-metadata-2.0.0-20240415.gdb"
    File gambit_db_signatures = "gs://gambit-databases-rp/2.0.0/gambit-signatures-2.0.0-20240415.gs"
    Int disk_size = 100
    Int memory = 2
    Int cpu = 1
  }
  command <<<
    # capture date
    date | tee DATE
    
    # copy database to a findable location within the execution directory
    # also, capture the version of the database used
    gambit_db_version="$(basename -- '~{gambit_db_genomes}'); $(basename -- '~{gambit_db_signatures}')"
    gambit_db_dir="${PWD}/gambit_database"
    mkdir ${gambit_db_dir}
    cp ~{gambit_db_genomes} ${gambit_db_dir}
    cp ~{gambit_db_signatures} ${gambit_db_dir}

    echo ${gambit_db_version} | tee GAMBIT_DB_VERSION
    
    # run gambit core check on the assembly
    echo "Running gambit core check with: gambit-core-check -e -v ${gambit_db_dir} ${gambit_db_dir}/$(basename -- '~{gambit_db_signatures}') ${gambit_db_dir}/$(basename -- '~{gambit_db_genomes}') ~{assembly}"
    gambit-core-check -e -v ${gambit_db_dir} ${gambit_db_dir}/$(basename -- '~{gambit_db_signatures}') ${gambit_db_dir}/$(basename -- '~{gambit_db_genomes}') ~{assembly} | tee ~{samplename}_gambit_core_report.tsv

    # parse output file
    cat ~{samplename}_gambit_core_report.tsv | cut -f 2 | tail -n 1 | tee SPECIES
    cat ~{samplename}_gambit_core_report.tsv | cut -f 3 | tail -n 1 | tee COMPLETNESS
    cat ~{samplename}_gambit_core_report.tsv | cut -f 4 | tail -n 1 | tee CORE_KMERS
    cat ~{samplename}_gambit_core_report.tsv | cut -f 5 | tail -n 1 | tee CLOSEST_ACCESSION
    cat ~{samplename}_gambit_core_report.tsv | cut -f 6 | tail -n 1 | tee CLOSEST_DISTANCE

  >>>
  output {
    File gambit_core_check_report_file = "~{samplename}_gambit_core_report.tsv"
    String gambit_core_check_species = read_string("SPECIES")
    String gambit_core_check_completeness = read_string("COMPLETNESS")
    String gambit_core_check_core_kmers = read_string("CORE_KMERS")
    String gambit_core_check_closest_accession = read_string("CLOSEST_ACCESSION")
    String gambit_core_check_closest_distance = read_string("CLOSEST_DISTANCE")
    String gambit_core_check_db_version = read_string("GAMBIT_DB_VERSION")
    String gambit_core_check_docker = docker
  }
  runtime {
    docker:  "~{docker}"
    memory:  "~{memory} GB"
    cpu:   "~{cpu}"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible:  0
  }
}