version 1.0

task metabuli {
  input {
    File read1 # intended for ONT reads only (at this time)
    String read1_basename = basename(basename(basename(read1, ".gz"), ".fastq"), ".fq")
    String samplename
    String taxon_of_interest

    #File metabuli_db = "gs://theiagen-large-public-files-rp/terra/databases/metabuli/refseq_virus-v223.tar.gz"
    File metabuli_db

    #File taxonomy_path = "gs://theiagen-large-public-files-rp/terra/databases/metabuli/new_taxdump.tar.gz"
    File taxonomy_path

    Float? min_score # metabuli: Min. sequence similarity score (0.0-1.0) [0.000]
    Float? min_sp_score # metabuli: Min. score for species- or lower-level classification. [0.000]
    Float? min_cov # metabuli: Min. query coverage (0.0-1.0) [0.000]
    Int cpu = 2
    Int memory = 4
    Int disk_size = 100
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/metabuli:1.1.0"
  }
  command <<<
    set -euo pipefail

    # Decompress additional taxonomy files necessary for ncbi ref_seq database search
    mkdir taxdump
    tar -C taxdump/ -xzvf ~{taxonomy_path}

    # Decompress/extract the ref_seq viral database
    mkdir db
    tar -C db/ -xzvf ~{metabuli_db}
    extracted_db=$(ls -d db/*/ | head -n 1)

    metabuli classify \
      --seq-mode 3 \
      ~{read1} \
      ${extracted_db} \
      output_dir \
      "~{samplename}" \
      --taxonomy-path taxdump/ \
      ~{"--min-score " + min_score} \
      ~{"--min-sp-score " + min_sp_score} \
      ~{"--min-cov " + min_cov} \
      --threads ~{cpu} \
      --max-ram ~{memory}

    metabuli extract \
      --seq-mode 3 \
      ~{read1} \
      ./output_dir/~{samplename}_classifications.tsv \
      ${extracted_db} \
      ~{"--tax-id " + taxon_of_interest}

    # the extracted reads are being output to the _miniwdl_input directory for some reason
    # my guess is metabuli is written in c++ and it's allocating memory for the output file before it's being executed
    # I don't think this will be needed when running on terra...
    find . -type f -name "*.fq" -exec mv {} . \;
  >>>
  output {
    File metabuli_report = "output_dir/~{samplename}_report.tsv"
    File metabuli_classified = "output_dir/~{samplename}_classifications.tsv"
    File metabuli_read1_extract = "~{read1_basename}_~{taxon_of_interest}.fq"
    String metabuli_docker = "~{docker_image}"
    String metabuli_database = "~{metabuli_db}"
  }
  runtime {
    docker: docker_image
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 0
  }
}