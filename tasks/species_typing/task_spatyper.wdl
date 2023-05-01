version 1.0

task spatyper {
  meta {
    description: "Computational method for finding spa types in Staphylococcus aureus"
  }
  input {
    File assembly
    String samplename
    String docker = "quay.io/biocontainers/spatyper:0.3.3--pyhdfd78af_3"
    Int disk_size = 100
    Int cpu = 4

    # Parameters
    # --do_enrich Do PCR product enrichment
    Boolean do_enrich = false
  }
  command <<<
    # get versioning
    spaTyper --version 2>&1 | sed 's/^.*spaTyper //' | tee VERSION
    
    spaTyper \
      ~{true="--do_enrich" false="" do_enrich} \
      --fasta ~{assembly} \
      --output ~{samplename}.tsv

    python3 <<CODE
    import csv

    TYPE = []
    REPEATS = []

    with open("./~{samplename}.tsv",'r') as tsv_file:
      tsv_reader=csv.reader(tsv_file, delimiter="\t")
      next(tsv_reader, None)  # skip the headers
      for row in tsv_reader:
        TYPE.append(row[-1])
        REPEATS.append(row[-2])

      with open ("TYPE", 'wt') as TYPE_fh:
        TYPE_fh.write(','.join(TYPE))

      with open ("REPEATS", 'wt') as REPEATS_fh:
        REPEATS_fh.write(','.join(REPEATS))
    CODE
  >>>
  output {
    File spatyper_tsv = "~{samplename}.tsv"
    String spatyper_repeats = read_string("REPEATS")
    String spatyper_type = read_string("TYPE")
    String spatyper_version = read_string("VERSION")
    String spatyper_docker = "~{docker}"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}
