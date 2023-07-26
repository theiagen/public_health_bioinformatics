version 1.0

task consensus {
  input {
    String samplename
    String? organism
    File filtered_reads
    File primer_bed
    File? reference_genome
    Int normalise = 20000
    Int cpu = 8
    Int disk_size = 100
    String medaka_model = "r941_min_high_g360"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/artic-ncov2019-epi2me"
  }
  String primer_name = basename(primer_bed)
  command <<<
    # HIV
    if [[ ~{organism} == "HIV" ]]; then
      # setup custom primer scheme (/V is required by Artic)
      mkdir -p ./primer-schemes/HIV/Vuser

      ## set reference genome
      ref_genome="~{reference_genome}"

      head -n1 "${ref_genome}" | sed 's/>//' | tee REFERENCE_GENOME
      cp "${ref_genome}" ./primer-schemes/HIV/Vuser/HIV.reference.fasta

      ## set primers
      #cp ~{primer_bed} ./primer-schemes/SARS-CoV-2/Vuser/SARS-CoV-2.scheme.bed
      #p_bed="~{primer_bed}"
      cp "~{primer_bed}" ./primer-schemes/HIV/Vuser/HIV.scheme.bed
      scheme_name="HIV/Vuser"
    # Add other viruses here
    # Default is SARS-CoV-2
    else
      # setup custom primer scheme (/V is required by Artic)
      mkdir -p ./primer-schemes/SARS-CoV-2/Vuser

      ## set reference genome
      if [[ ! -z "~{reference_genome}" ]]; then
        ref_genome="~{reference_genome}"
      else
        # use reference file in docker--different paths depending on image specified 
        if [[ -d "/fieldbioinformatics" ]]; then
          ref_genome=$(find /fieldbioinformatics/*/primer*schemes/nCoV-2019/V3/ -name "nCoV-2019.reference.fasta")
        else
          ref_genome=$(find /wf-artic*/data/primer_schemes/SARS-CoV-2/V4/ -name "SARS-CoV-2.reference.fasta")
        fi
        echo "No user-defined reference genome; setting reference to ${ref_genome}"
      fi
      head -n1 "${ref_genome}" | sed 's/>//' | tee REFERENCE_GENOME
      cp "${ref_genome}" ./primer-schemes/SARS-CoV-2/Vuser/SARS-CoV-2.reference.fasta

      ## set primers
      cp ~{primer_bed} ./primer-schemes/SARS-CoV-2/Vuser/SARS-CoV-2.scheme.bed
      scheme_name="SARS-CoV-2/Vuser"
    fi

    # version control
    echo "Medaka via $(artic -v)" | tee VERSION
    echo "~{primer_name}" | tee PRIMER_NAME
    artic minion --medaka --medaka-model ~{medaka_model} --normalise ~{normalise} --threads ~{cpu} --scheme-directory ./primer-schemes --read-file ~{filtered_reads} ${scheme_name} ~{samplename}
    gunzip -f ~{samplename}.pass.vcf.gz

    # clean up fasta header
    echo ">~{samplename}" > ~{samplename}.medaka.consensus.fasta
    grep -v ">" ~{samplename}.consensus.fasta >> ~{samplename}.medaka.consensus.fasta

    # grab reads from alignment
    samtools fastq -F4 ~{samplename}.primertrimmed.rg.sorted.bam | gzip > ~{samplename}.fastq.gz  
  >>>
  output {
    File consensus_seq = "~{samplename}.medaka.consensus.fasta"
    File sorted_bam = "~{samplename}.trimmed.rg.sorted.bam"
    File trim_sorted_bam = "~{samplename}.primertrimmed.rg.sorted.bam"
    File trim_sorted_bai = "~{samplename}.primertrimmed.rg.sorted.bam.bai"
    File medaka_pass_vcf = "~{samplename}.pass.vcf"
    File? reads_aligned = "~{samplename}.fastq.gz"
    String medaka_reference = read_string("REFERENCE_GENOME")
    String artic_pipeline_version = read_string("VERSION")
    String artic_pipeline_docker = docker
    String primer_bed_name = read_string("PRIMER_NAME")
    File? trim_fastq = "~{samplename}.primertrimmed.rg.fastq"
  }
  runtime {
    docker: "~{docker}"
    memory: "16 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
    maxRetries: 3
  }
}