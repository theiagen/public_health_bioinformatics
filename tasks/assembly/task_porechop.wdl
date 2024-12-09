version 1.0

task porechop {
  input {
    Array[File] raw_reads             # Input raw reads (FASTQ format)
    String samplename                 
    Int cpu = 4                     
    Int memory = 8                    
    Int disk_size = 50   
    String docker #need to create dockerfile      
  }

  command <<< 
        porechop_abi --version | tee VERSION

        # Combine input reads into a single FASTQ file (if more than one file)
        cat ~{sep=" " raw_reads} > combined_raw_reads.fastq

        # Run Porechop
        porechop_abi \
          -input combined_raw_reads.fastq \
          -output ~{samplename}.trimmed.fastq \
          --threads ~{cpu} 
    >>>

  output {
    File trimmed_reads = "~{samplename}.trimmed.fastq"   # Cleaned reads
    String porechopabi_version = read_string("VERSION")    # Porechop version
  }
  runtime {
    docker: "~{docker}"
    cpu: cpu
    memory: "~{memory} GB"
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }   
}
