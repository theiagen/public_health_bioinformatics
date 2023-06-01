version 1.0

task compare_assemblies {
  meta {
    description: "Compares the de novo and consensus assembly and returns the one with the highest base count"
  }
  input {
    File assembly_denovo
    File assembly_consensus
    String samplename
    Int disk_size = 100
  }
  command <<<
    python3 <<CODE
    import os

    # Auxiliary functions
    def file_is_empty(path):
        return os.stat(path).st_size==0

    def count_bases(file_name):
        """Count the number of bases in a FASTA file."""
        with open(file_name, "r") as f:
            seq = ""
            for line in f:
                if line.startswith(">"):
                    continue
                seq += line.strip()
            seq = seq.replace("N","")  # ignore uncalled bases
            return len(seq)
    
    # File comparison
    denovo_base_count = count_bases("~{assembly_denovo}")
    consensus_base_count = count_bases("~{assembly_consensus}")


    if denovo_base_count >= consensus_base_count:
        input_file = "~{assembly_denovo}"
    else:
        input_file = "~{assembly_consensus}"

    with open("~{samplename}_highest.fasta", 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                out_f.write(line)
    
    # File validations
    if file_is_empty("~{samplename}_highest.fasta"):
        print("Assembly file is empty! Removing it...")
        os.remove("~{samplename}_highest.fasta")
    
    if input_file == "~{assembly_consensus}" and consensus_base_count == 0:
        print("Assembly file contains just uncalled bases! Removing it...")
        os.remove("~{samplename}_highest.fasta")

    CODE
  >>>
  output {
    File final_assembly = "~{samplename}_highest.fasta"
  }
  runtime {
    docker: "quay.io/theiagen/terra-tools:2023-03-16"
    memory: "8 GB"
    cpu: 1
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 0
  }
}