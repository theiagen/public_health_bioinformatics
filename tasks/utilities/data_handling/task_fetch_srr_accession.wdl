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
    version: "1.0"
  }

  command <<< 
    date -u | tee DATE

    fastq-dl --version | tee VERSION

    # Debug output to show the sample being processed
    echo "Fetching metadata for sample accession: ${sample_accession}"

    # Use fastq-dl to fetch metadata only, outputting to the current directory
    fastq-dl --accession ~{sample_accession} --only-download-metadata --verbose

    if [[ -f fastq-run-info.tsv ]]; then
      echo "Metadata written for ${sample_accession}:"
      echo "TSV content:"
      cat fastq-run-info.tsv

      # Extract the SRR accessions and write them directly to srr_accession.txt
      awk -F'\t' 'NR>1 {print $1}' fastq-run-info.tsv > srr_accession.txt

      # Check if srr_accession.txt is empty (no SRR accessions found)
      if [[ ! -s srr_accession.txt ]]; then
        echo "No SRR accession found for ${sample_accession}" > srr_accession.txt
      else
        echo "Extracted SRR accessions:"
        cat srr_accession.txt
      fi
    else
      echo "No metadata found for ${sample_accession}"
      echo "No SRR accession found" > srr_accession.txt
    fi
  >>>
  output {
    Array[String] srr_accession = read_lines("srr_accession.txt")
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
