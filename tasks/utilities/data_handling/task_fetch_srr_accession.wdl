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

    # Fetch metadata for the sample accession
    echo "Fetching metadata for valid biosample ID or SRA: ~{sample_accession}"
    if fastq-dl --accession ~{sample_accession} --only-download-metadata --verbose; then
        if [[ -f fastq-run-info.tsv ]]; then
            echo "Metadata written for valid biosample ID or SRA: ~{sample_accession}"
            cat fastq-run-info.tsv

            # Extract SRR accessions from the TSV file
            SRR_accessions=$(awk -F'\t' 'NR>1 {print $1}' fastq-run-info.tsv | paste -sd ',' -)

            if [[ -z "${SRR_accessions}" ]]; then
                echo "No SRR accession found for valid biosample ID or SRA: ~{sample_accession}" > srr_accession.txt
            else
                echo "Extracted SRR accessions: ${SRR_accessions}"
                echo "${SRR_accessions}" > srr_accession.txt
            fi
        else
            echo "No metadata file found for valid biosample ID or SRA: ~{sample_accession}"
            echo "No SRR accession found" > srr_accession.txt
        fi
    else
        # Handle cases where fastq-dl exits with an error
        if grep -q "No results found" stderr || [[ ! -f fastq-run-info.tsv ]]; then
            echo "No SRR accession found for valid biosample ID or SRA: ~{sample_accession}"
            echo "No SRR accession found" > srr_accession.txt
        else
            echo "fastq-dl failed for ~{sample_accession}, invalid biosample ID or SRA"
            exit 1
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
