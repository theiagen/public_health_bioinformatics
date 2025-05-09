version 1.0

task extract_kraken_reads {
  input {
    File kraken2_output
    File kraken2_report
    File read1
    File read2
    Int taxon_id

    Int cpu = 1
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/krakentools:d4a2fbe"
    Int memory = 4
  }
  command <<<
    # fail hard
    set -euo pipefail

    # decompress classified data if it is gzipped
    gunzip -c ~{kraken2_output} > kraken2_output_unzipped.txt

    # add exclusion tag
    echo "DEBUG: Extracting reads"
    python3 /KrakenTools/extract_kraken_reads.py \
      -k kraken2_output_unzipped.txt \
      -s1 ~{read1} \
      -s2 ~{read2} \
      --taxid ~{taxon_id} \
      --report ~{kraken2_report} \
      --include-children \
      --fastq-output \
      --output ~{taxon_id}_1.fastq \
      --output2 ~{taxon_id}_2.fastq \

    if [ -s ~{taxon_id}_1.fastq ]; then
      echo "DEBUG: Taxon ~{taxon_id} reads extracted"
      echo "true" > CONTINUE

      gzip ~{taxon_id}_1.fastq 
      gzip ~{taxon_id}_2.fastq
    else
      echo "DEBUG: No reads were extracted for taxon ~{taxon_id}, removing empty files"
      echo "false" > CONTINUE
    fi
    
    grep ~{taxon_id} ~{kraken2_report} | awk '{for (i=6; i <= NF; ++i) print $i}' | tr '\n' ' ' | xargs > ORGANISM_NAME
  >>>
  output {
    File? extracted_read1 = "~{taxon_id}_1.fastq.gz"
    File? extracted_read2 = "~{taxon_id}_2.fastq.gz"
    String organism_name = read_string("ORGANISM_NAME")
    String krakentools_docker = docker
    Boolean success = read_boolean("CONTINUE")
  }
  runtime {
    cpu: cpu
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: docker
    memory: "~{memory} GB"
    preemptible: 1
    maxRetries: 3
  }
}