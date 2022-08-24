version 1.0

task mycosnp {
  input {
    File read1
    File read2
    String samplename
    String docker = "quay.io/theiagen/mycosnp:dev"
    String strain = "B11205"
    String accession = "GCA_016772135"
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    date | tee DATE
    echo $(nextflow pull rpetit3/mycosnp-nf 2>&1) | sed 's/^.*revision: //;' | tee MYCOSNP_VERSION

    # Make sample FOFN
    echo "sample,fastq_1,fastq_2" > sample.csv
    echo "~{samplename},~{read1},~{read2}" >> sample.csv

    # Run MycoSNP
    mkdir ~{samplename}
    cd ~{samplename}
    if nextflow run rpetit3/mycosnp-nf --input ../sample.csv --ref_dir /reference/~{accession} --publish_dir_mode copy --skip_phylogeny; then
      # Everything finished, pack up the results and clean up
      rm -rf .nextflow/ work/
      cd ..
      tar -cf - ~{samplename}/ | gzip -n --best  > ${samplename}.tar.gz
    else
      # Run failed
      exit 1
    fi
  >>>
  output {
    String mycosnp_version = read_string("MYCOSNP_VERSION")
    String mycosnp_docker = docker
    String analysis_date = read_string("DATE")
    String reference_strain = strain
    String reference_accession = accession
    File assembly_fasta = "~{samplename}/results/combined/consensus/~{samplename}.fasta.gz"
    File full_results = "~{samplename}.tar.gz"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk 50 SSD"
    maxRetries: 3
    preemptible: 0
  }
}
