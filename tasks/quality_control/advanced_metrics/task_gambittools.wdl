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
    
    gambit_db_version="$(basename -- '~{gambit_db_genomes}'); $(basename -- '~{gambit_db_signatures}')"
    gambit_db_dir="${PWD}/gambit_database"
    mkdir ${gambit_db_dir}
    cp ~{gambit_db_genomes} ${gambit_db_dir}
    cp ~{gambit_db_signatures} ${gambit_db_dir}

    echo ${gambit_db_version} | tee GAMBIT_DB_VERSION
    
    echo "Running gambit core check with: gambit-core-check -e -v ${gambit_db_dir} ${gambit_db_dir}/$(basename -- '~{gambit_db_signatures}') ${gambit_db_dir}/$(basename -- '~{gambit_db_genomes}') ~{assembly}"
    gambit-core-check -e -v ${gambit_db_dir} ${gambit_db_dir}/$(basename -- '~{gambit_db_signatures}') ${gambit_db_dir}/$(basename -- '~{gambit_db_genomes}') ~{assembly} | tee ~{samplename}_gambit_core_report.txt
  >>>
  output {
    File gambit_core_check_report_file = "~{samplename}_gambit_core_report.txt"
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