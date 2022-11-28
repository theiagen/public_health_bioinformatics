version 1.0

workflow ecoli_char {

  input {
    String    SRR
    File      contigs
  }

  call abricate as abricate {
    input:
      samplename=SRR,
      contigs=contigs,
      database="ncbi"
  }

  call abricate as abricate_virfinder {
    input:
      samplename=SRR,
      contigs=contigs,
      database="ecoli_vf"
  }

  call amrfinderplus {
    input:
      samplename=SRR,
      contigs=contigs
  }

  call serotypefinder {
    input:
      samplename=SRR,
      contigs=contigs
  }

  output {
    File    abricate_results                =abricate.abricate_results
    File    abricate_virfinder_results                =abricate_virfinder.abricate_results
    File    amrfinderplus_results           =amrfinderplus.amrfinder_results
    File    serotypefinder_results          =serotypefinder.serotypefinder_results
  }
}

task abricate {

  input {
    File      contigs
    String    samplename
    String    database
  }

  command {
    abricate --version | head -1 | tee VERSION
    abricate --db ${database} ${contigs} > ${samplename + '_abricate.tsv'}
  }

  output {
    File      abricate_results="${samplename + '_abricate.tsv'}"
  }

  runtime {
    docker:       "quay.io/staphb/abricate:1.0.0"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task amrfinderplus {
  input {
    File      contigs
    String    samplename
  }

  command {
    amrfinder --version | head -1 | tee VERSION
    amrfinder \
    --nucleotide ${contigs} \
    -o ${samplename + '_amrfinder.tsv'}
  }

  output {
    File      amrfinder_results="${samplename + '_amrfinder.tsv'}"
  }

  runtime {
    docker:       "quay.io/staphb/ncbi-amrfinderplus:3.8.28"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task serotypefinder {

  input {
    File      contigs
    String    samplename
  }

  command {
    serotypefinder.pl --version | head -1 | tee VERSION
    serotypefinder.pl \
    -i ${contigs} \
    -d /serotypefinder/database \
    -b /blast-2.2.26 \
    -s ecoli \
    -k 85.00 \
    -l 0.60 \
    -o ${samplename}
  }

  output {
    File      serotypefinder_results="${samplename}/results_table.txt"
  }

  runtime {
    docker:       "quay.io/staphb/serotypefinder:1.1"
    memory:       "8 GB"
    cpu:          2
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
