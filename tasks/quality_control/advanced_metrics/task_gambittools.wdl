version 1.0

task gambitcore {
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/internal/gambitcore:0.0.3"
    File gambit_db_genomes = "gs://gambit-databases-rp/1.3.0/gambit-metadata-1.3-231016.gdb"
    File gambit_db_signatures = "gs://gambit-databases-rp/1.3.0/gambit-signatures-1.3-231016.gs"
    Int disk_size = 100
    Int memory = 2
    Int cpu = 1
    Int? max_species_genomes # Max number of genomes in a species to consider, ignore all others above this
    Float? core_proportion # Proportion of genomes a kmer must be in for a species to be considered core
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
    echo "Running gambitcore with: gambitcore ${gambit_db_dir} ~{assembly}"
    gambitcore \
      --cpus ~{cpu} \
      ~{'--max_species_genomes ' + max_species_genomes} \
      ~{'--core_proportion ' + core_proportion} \
      ${gambit_db_dir} ~{assembly} | tee ~{samplename}_gambitcore_report.tsv
    
    # check if the output file was created
    if [ ! -s ~{samplename}_gambitcore_report.tsv ]; then
      echo "Gambitcore failed to produce a report"
      exit 1
    fi

    # check if file has more than 1 line
    if [ $(wc -l < ~{samplename}_gambitcore_report.tsv) -lt 2 ]; then
      echo "Gambitcore produced an empty report"
      exit 1
    fi

    # parse output file
    cat ~{samplename}_gambitcore_report.tsv | cut -f 2 | tail -n 1 | tee SPECIES
    cat ~{samplename}_gambitcore_report.tsv | cut -f 3 | tail -n 1 | tee COMPLETENESS
    cat ~{samplename}_gambitcore_report.tsv | cut -f 4 | tail -n 1 | tee N_GENOMES
    cat ~{samplename}_gambitcore_report.tsv | cut -f 5 | tail -n 1 | tee ASSEMBLY_SPECIES_CORE_KMERS
    cat ~{samplename}_gambitcore_report.tsv | cut -f 6 | tail -n 1 | tee CLOSEST_ACCESSION
    cat ~{samplename}_gambitcore_report.tsv | cut -f 7 | tail -n 1 | tee CLOSEST_DISTANCE
    cat ~{samplename}_gambitcore_report.tsv | cut -f 8 | tail -n 1 | tee ASSEMBLY_KMERS
    cat ~{samplename}_gambitcore_report.tsv | cut -f 9 | tail -n 1 | tee SPECIES_KMERS
    cat ~{samplename}_gambitcore_report.tsv | cut -f 10 | tail -n 1 | tee SPECIES_STD_KMERS
    cat ~{samplename}_gambitcore_report.tsv | cut -f 11 | tail -n 1 | tee ASSEMBLY_QC

  >>>
  output {
    File gambitcore_report_file = "~{samplename}_gambitcore_report.tsv"
    String gambitcore_species = read_string("SPECIES")
    String gambitcore_completeness = read_string("COMPLETENESS")
    String gambitcore_n_genomes = read_string("N_GENOMES")
    String gambitcore_kmers_ratio = read_string("ASSEMBLY_SPECIES_CORE_KMERS")
    String gambitcore_closest_accession = read_string("CLOSEST_ACCESSION")
    String gambitcore_closest_distance = read_string("CLOSEST_DISTANCE")
    String gambitcore_assembly_kmers = read_string("ASSEMBLY_KMERS")
    String gambitcore_species_kmers = read_string("SPECIES_KMERS")
    String gambitcore_species_std_kmers = read_string("SPECIES_STD_KMERS")
    String gambitcore_assembly_qc = read_string("ASSEMBLY_QC")
    String gambitcore_db_version = read_string("GAMBIT_DB_VERSION")
    String gambitcore_docker = docker
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