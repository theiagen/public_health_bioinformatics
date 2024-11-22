version 1.0

task dnaapler_all {
  input {
    File input_fasta                       # Input genome assembly in FASTA format
    String output_dir                      # Directory to store output
    String prefix       
    Int threads = 4         
    Int disk_size = 100               
    Int memory = 16                        
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/dnaapler:0.8.0"
 }

  command <<< 
      dnaapler --version | tee VERSION

      # Run DnaApler with the 'all' subcommand
      dnaapler all \
        -i ~{input_fasta} \
        -o ~{output_dir} \
        -p ~{prefix} \
        -t ~{threads}
  >>>

  output {
    File reoriented_fasta = "~{output_dir}/~{prefix}_reoriented.fasta"  # Reoriented genome assembly
    String dnaapler_version = read_string("VERSION")                   # DnaApler version
  }
  runtime {
    docker: "~{docker}"           
    cpu: threads                       
    memory: "~{memory} GB"                  
    disks: "local-disk " + disk_size + " SSD"  
    disk: disk_size + " GB"             
    maxRetries: 3                     
    preemptible: 0                     
  }

