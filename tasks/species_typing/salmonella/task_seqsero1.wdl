version 1.0

task seqsero1_assembly {
  input {
    File assembly_fasta
    String samplename
    String docker = "staphb/seqsero:1.0.1"
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    set -euo pipefail

    # Print and save version
    echo "SeqSero v1.0.1" | tee VERSION

    SeqSero.py \
      -m 4 \
      -d ~{samplename}_seqsero_output_dir \
      -i ~{assembly_fasta}

    touch PREDICTED_ANTIGENIC_PROFILE PREDICTED_SEROTYPE

    grep "Predicted antigenic profile:" ~{samplename}_seqsero_output_dir/Seqsero_result.txt | cut -f2- | tee PREDICTED_ANTIGENIC_PROFILE
    grep "Predicted serotype(s):" ~{samplename}_seqsero_output_dir/Seqsero_result.txt | cut -f2- | tee PREDICTED_SEROTYPE

  >>>
  output {
    File seqsero_report = "./~{samplename}_seqsero_output_dir/Seqsero_result.txt"
    String seqsero_version = read_string("VERSION")
    String seqsero_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero_predicted_serotype = read_string("PREDICTED_SEROTYPE")
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}