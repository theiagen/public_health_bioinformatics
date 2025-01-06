version 1.0

task consensus {
  input {
    String samplename
    String? organism
    File read1
    File? primer_bed
    File? reference_genome
    Int normalise = 20000
    Int cpu = 8
    Int memory = 16
    Int disk_size = 100
    String clair3_model = "r941_prom_hac_g360+g422"
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/artic:1.6.0"
  }
  command <<<
    set -euo pipefail

    # Determine if we're using provided or local schemes - most of the time will be provided
    # But with newer versions of ARTIC, we can use remote schemes from https://github.com/quick-lab/primerschemes
    if [[ -z "~{primer_bed}" && -z "~{reference_genome}" ]]; then
      echo "Using remote scheme..."
      
      # Set scheme parameters based on organism
      if [[ "~{organism}" == "sars-cov-2" ]]; then
        scheme_name="sars-cov-2"
        scheme_version="v4.0.0"
      elif [[ "~{organism}" == "MPXV" ]]; then
        scheme_name="mpox"
        scheme_version="v1.0.0"
      else
        echo "Error: Unsupported organism for remote schemes: ~{organism}" >&2
        echo "Please provide primer bed and reference genome for custom organisms" >&2
        exit 1
      fi

      # Run ARTIC with remote scheme per newer versions of ARTIC
      artic minion --model ~{clair3_model} \
        --normalise ~{normalise} \
        --threads ~{cpu} \
        --scheme-directory /data/primer-schemes \
        --scheme-name ${scheme_name} \
        --scheme-version ${scheme_version} \
        --read-file ~{read1} \
        ~{samplename}

      # Set output file contents for remote scheme
      echo "${scheme_name}" > REFERENCE_GENOME
      echo "${scheme_name}_${scheme_version}" > PRIMER_NAME

    else
      # Newer version of Artic use a different command to launch the pipeline
      # when using user provided files for bed and reference genome are provided, 
      # we can now provide the bed and reference genome files to the pipeline directly 
      # instead of using the --scheme-directory and --scheme-name flags
      echo "Using user provided files..."
      
      if [[ -z "~{primer_bed}" || -z "~{reference_genome}" ]]; then
        echo "Error: Both primer bed and reference genome must be provided for local scheme" >&2
        exit 1
      fi

      # Run ARTIC with user provided files
      artic minion --model ~{clair3_model} \
        --normalise ~{normalise} \
        --threads ~{cpu} \
        --bed ~{primer_bed} \
        --ref ~{reference_genome} \
        --read-file ~{read1} \
        ~{samplename}

      # Set output file contents for local scheme
      head -n1 "~{reference_genome}" | sed 's/>//' > REFERENCE_GENOME
      basename "~{primer_bed}" > PRIMER_NAME
    fi

    # Capture ARTIC version
    echo "Artic Pipeline Version $(artic -v)" > VERSION

    # Grab reads from alignment - cdph wants this
    samtools fastq -F4 ~{samplename}.primertrimmed.rg.sorted.bam | gzip > ~{samplename}.fastq.gz  
  >>>
  output {
    File artic_consensus_fasta = "~{samplename}.consensus.fasta"
    File artic_clair3_pass_vcf = "~{samplename}.pass.vcf"
    String artic_pipeline_reference = read_string("REFERENCE_GENOME")
    String artic_pipeline_version = read_string("VERSION")
    String artic_pipeline_docker = docker
    File sorted_bam = "~{samplename}.trimmed.rg.sorted.bam"
    File trim_sorted_bam = "~{samplename}.primertrimmed.rg.sorted.bam"
    File trim_sorted_bai = "~{samplename}.primertrimmed.rg.sorted.bam.bai"
    File? reads_aligned = "~{samplename}.fastq.gz"
    String primer_bed_name = read_string("PRIMER_NAME")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}