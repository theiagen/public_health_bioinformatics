version 1.0

task fetch_srr_accession {
  input {
    String sample_accession 
    String docker = "us-docker.pkg.dev/general-theiagen/biocontainers/fastq-dl:2.0.4--pyhdfd78af_0"
    Int disk_size = 10
    Int cpu = 2
    Int memory = 8
  }
  meta {
    volatile: true
  }
  command <<< 
    set -euo pipefail

    # Output the current date and fastq-dl version for debugging
    date -u | tee DATE
    fastq-dl --version | tee VERSION

    echo "Fetching metadata for accession: ~{sample_accession}"

    # Run fastq-dl and capture stderr
    fastq-dl --accession ~{sample_accession} --only-download-metadata -m 2 --verbose 2> stderr.log || true

    # Handle whether the ID/accession is valid and contains SRR metadata based on stderr
    if grep -q "No results found for" stderr.log; then
        echo "No SRR accession found" > srr_accession.txt
        echo "No SRR accession found for accession: ~{sample_accession}"
    elif grep -q "received an empty response" stderr.log; then
        echo "No SRR accession found" > srr_accession.txt
        echo "No SRR accession found for accession: ~{sample_accession}"
    elif grep -q "is not a Study, Sample, Experiment, or Run accession" stderr.log; then
        echo "Invalid accession: ~{sample_accession}" >&2
        exit 1
    elif [[ ! -f fastq-run-info.tsv ]]; then
        echo "No metadata file found for accession: ~{sample_accession}" >&2
        exit 1
    else
        # Extract SRR accessions from the TSV file if it exists
        SRR_accessions=$(awk -F'\t' 'NR>1 {print $1}' fastq-run-info.tsv | paste -sd ',' -)
        if [[ -z "${SRR_accessions}" ]]; then
            echo "No SRR accession found" > srr_accession.txt
        else
            echo "Extracted SRR accessions: ${SRR_accessions}"
            echo "${SRR_accessions}" > srr_accession.txt
        fi
    fi
  >>>
  output {
    String srr_accession = read_string("srr_accession.txt")
    String fastq_dl_version = read_string("VERSION")
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}
