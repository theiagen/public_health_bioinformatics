version 1.0

task sonneityping {
  # Inputs
  input {
    File read1
    File? read2
    Boolean ont_data = false
    String samplename
    String docker = "staphb/mykrobe:0.12.1"
    Int disk_size = 100
    String? myrkobe_opts
    Int cpu = 4
  }
  command <<<
    # Print and save versions
     mykrobe --version | sed 's|mykrobe v||g' | tee MYKROBE_VERSION.txt
    # opting to skip capturing the sonneityping version since there is no --version flag or easy way to determine version
    # navigate here for docker image and version information: https://github.com/StaPH-B/docker-builds/tree/master/mykrobe

    # Run Mykrobe on the input read data
    mykrobe predict \
    -t ~{cpu} \
    --sample ~{samplename} \
    --species sonnei \
    --format json_and_csv \
    --out ~{samplename}.mykrobe \
    ~{true='--ont' false='' ont_data} \
    --seq ~{read1} ~{read2} \
    ~{myrkobe_opts}

    # use sonneityping script to produce final TSV; alleles.txt is required input for human-readable genotype names
    python /sonneityping/parse_mykrobe_predict.py \
    --jsons ~{samplename}.mykrobe.json --alleles /sonneityping/alleles.txt \
    --prefix ~{samplename}.sonneityping

    # rename output TSV to something prettier
    mv -v ~{samplename}.sonneityping_predictResults.tsv ~{samplename}.sonneityping.tsv

    # Run a python block to parse output sonneityping TSV file for terra data tables
    python3 <<CODE
    import csv
    with open("./~{samplename}.sonneityping.tsv",'r') as tsv_file:
      tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
      for line in tsv_reader:
        with open ("SPECIES.txt", 'wt') as sonneityping_species:
          species=line["species"]
          sonneityping_species.write(species)
        with open ("FINAL_GENOTYPE.txt", 'wt') as final_genotype:
          genotype=line["final genotype"]
          final_genotype.write(genotype)
        with open ("GENOTYPE_NAME.txt", 'wt') as genotype_name:
          genotypename=line["name"]
          genotype_name.write(genotypename)
        with open ("CONFIDENCE.txt", 'wt') as sonneityping_confidence:
          confidence=line["confidence"]
          sonneityping_confidence.write(confidence)
    CODE
  >>>
  output {
    File sonneityping_mykrobe_report_csv = "~{samplename}.mykrobe.csv"
    File sonneityping_mykrobe_report_json = "~{samplename}.mykrobe.json"
    File sonneityping_final_report_tsv = "~{samplename}.sonneityping.tsv"
    String sonneityping_mykrobe_version = read_string("MYKROBE_VERSION.txt")
    String sonneityping_mykrobe_docker = docker
    String sonneityping_species = read_string("SPECIES.txt")
    String sonneityping_final_genotype = read_string("FINAL_GENOTYPE.txt")
    String sonneityping_genotype_confidence = read_string("CONFIDENCE.txt")
    String sonneityping_genotype_name = read_string("GENOTYPE_NAME.txt")
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