version 1.0

task seqsero2 {
  # Inputs
  input {
    File read1
    File? read2
    String samplename
    String mode = "a"
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqsero2:1.3.1"
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
    Boolean paired_end
  }
  command <<<
    set -euo pipefail
    
    # Print and save version
    SeqSero2_package.py --version | tee VERSION

    # Run SeqSero2 on the input read data
    # explanation of parameters:
    # -p <int> (number of threads for allele mode)
    # -m <string> (which workflow to apply, either 'a' for raw reads allele micro-assembly 
    #              or 'k' for raw reads and genome assembly k-mer; default=a)
    # -t <int> (input data type; '2' for separated paired-end reads, '3' for single reads)
    # -n <string> (optional, to specify sample name in report output)
    # -d <string> (output directory name)
    # -i <file> (path to input file(s))
    SeqSero2_package.py \
      -p ~{cpu} \
      -m ~{mode} \
      ~{true='-t 2' false='-t 3' paired_end} \
      -n ~{samplename} \
      -d ~{samplename}_seqsero2_output_dir \
      -i ~{read1} ~{read2}

    # Run a python block to parse output file for terra data tables
    python3 <<CODE
    import csv
    with open("./~{samplename}_seqsero2_output_dir/SeqSero_result.tsv",'r') as tsv_file:
      tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
      
      for line in tsv_reader:
        with open ("PREDICTED_ANTIGENIC_PROFILE", 'wt') as Predicted_Antigen_Prof:
          pred_ant_prof=line['Predicted antigenic profile']
          if not pred_ant_prof:
            pred_ant_prof = "None"
          Predicted_Antigen_Prof.write(pred_ant_prof)
        
        with open ("PREDICTED_SEROTYPE", 'wt') as Predicted_Sero:
          pred_sero=line['Predicted serotype']
          if not pred_sero:
            pred_sero = "None"
          Predicted_Sero.write(pred_sero)
        
        with open ("CONTAMINATION", 'wt') as Contamination_Detected:
          cont_detect=line['Potential inter-serotype contamination']
          if not cont_detect:
            cont_detect = "None"
          Contamination_Detected.write(cont_detect)

        with open ("NOTE", 'wt') as Note:
          note=line['Note']
          if not note:
            note = "None"
          Note.write(note)
    CODE
  >>>
  output {
    File seqsero2_report = "./~{samplename}_seqsero2_output_dir/SeqSero_result.tsv"
    String seqsero2_version = read_string("VERSION")
    String seqsero2_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero2_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String seqsero2_predicted_contamination = read_string("CONTAMINATION")
    String seqsero2_note = read_string("NOTE")
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

task seqsero2_assembly {
  input {
    File assembly_fasta
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqsero2:1.3.1"
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    set -euo pipefail

    # Print and save version
    SeqSero2_package.py --version | tee VERSION

    # Run SeqSero2 on the input assembly fastas
    # explanation of parameters:
    # -p <int> (number of threads; 4 is the max that can be used for assembly data)
    # -m <string> ('k' for raw reads and genome assembly k-mer; required when t=4)
    # -t <int> (input data type; '4' for genome assembly)
    # -n <string> (optional, to specify sample name in report output)
    # -d <string> (output directory name)
    # -i <file> (path to input file)
    SeqSero2_package.py \
      -p 4 \
      -m k \
      -t 4 \
      -n ~{samplename} \
      -d ~{samplename}_seqsero2_output_dir \
      -i ~{assembly_fasta}

    # Run a python block to parse output file for terra data tables
    # contamination is not predicted when assembly data is used
    python3 <<CODE
    import csv
    with open("./~{samplename}_seqsero2_output_dir/SeqSero_result.tsv",'r') as tsv_file:
      tsv_reader = list(csv.DictReader(tsv_file, delimiter="\t"))
      
      for line in tsv_reader:
        with open ("PREDICTED_ANTIGENIC_PROFILE", 'wt') as Predicted_Antigen_Prof:
          pred_ant_prof=line['Predicted antigenic profile']
          if not pred_ant_prof:
            pred_ant_prof = "None"
          Predicted_Antigen_Prof.write(pred_ant_prof)
        
        with open ("PREDICTED_SEROTYPE", 'wt') as Predicted_Sero:
          pred_sero=line['Predicted serotype']
          if not pred_sero:
            pred_sero = "None"
          Predicted_Sero.write(pred_sero)
          
        with open ("NOTE", 'wt') as Note:
          note=line['Note']
          if not note:
            note = "None"
          Note.write(note)

    CODE
  >>>
  output {
    File seqsero2_report = "./~{samplename}_seqsero2_output_dir/SeqSero_result.tsv"
    String seqsero2_version = read_string("VERSION")
    String seqsero2_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero2_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String seqsero2_note = read_string("NOTE")
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