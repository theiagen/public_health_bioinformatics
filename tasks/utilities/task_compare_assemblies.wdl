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

    def count_bases(file_name):
        """Count the number of bases in a FASTA file."""
        with open(file_name, "r") as f:
            seq = ""
            for line in f:
                if line.startswith(">"):
                    continue
                seq += line.strip()
            return len(seq)
    
    denovo_base_count = count_bases("~{assembly_denovo}")
    consensus_base_count = count_bases("~{assembly_consensus}")

    if denovo_base_count > consensus_base_count:
        input_file = "~{assembly_denovo}"
    else:
        input_file = "~{assembly_consensus}"

    with open("~{samplename}_highest.fasta", 'w') as out_f:
        with open(input_file, 'r') as in_f:
            for line in in_f:
                out_f.write(line)

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
    maxRetries: 0 #TODO change to 3
  }
}