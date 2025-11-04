version 1.0

task seqsero2s {
  input {
    File read1
    File? read2
    String samplename
    String mode = "a"
    
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqsero2s:1.1.4"
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
    
    Boolean paired_end
  }
  command <<<
    set -euo pipefail
    
    # Print and save version
    SeqSero2S.py --version | tee VERSION

    # Usage: SeqSero2S.py 
    # -m <string> (which workflow to apply, 'a' (raw reads allele micro-assembly), 'k' (raw reads and genome assembly k-mer), default=a)
    # -t <string> (input data type; HARDCODED TO EITHER '2' (separated paired-end reads) OR '3' (single reads) SINCE ASSEMBLY INPUT IS HANDLED SEPARATELY)
    # -i <file> (/path/to/input/file)
    # -p <int> (number of threads for allele mode) 
    # -d <string> (output directory name)
    # -n <string> (sets the sample name in the report output)

    SeqSero2S.py \
      -m ~{mode} \
      -p ~{cpu} \
      ~{true='-t 2' false='-t 3' paired_end} \
      -n ~{samplename} \
      -d ~{samplename}_seqsero2s_output_dir \
      -i ~{read1} ~{read2}

    # Run a python block to parse output file for terra data tables
    python3 <<CODE
    import csv
    with open("./~{samplename}_seqsero2s_output_dir/SeqSero_result.tsv",'r') as tsv_file:
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
    File seqsero2s_report = "./~{samplename}_seqsero2s_output_dir/SeqSero_result.tsv"
    String seqsero2s_version = read_string("VERSION")
    String seqsero2s_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero2s_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String seqsero2s_predicted_contamination = read_string("CONTAMINATION")
    String seqsero2s_note = read_string("NOTE")
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

task seqsero2s_assembly {
  input {
    File assembly_fasta
    String samplename

    String docker = "us-docker.pkg.dev/general-theiagen/staphb/seqsero2s:1.1.4"
    Int disk_size = 100
    Int memory = 16
    Int cpu = 4
  }
  command <<<
    set -euo pipefail

    # Print and save version
    SeqSero2S.py --version | tee VERSION

    # Usage: SeqSero2S.py 
    # -m <string> (which workflow to apply; HARDCODED TO 'k' (raw reads and genome assembly k-mer) SINCE THE OTHER MODE DOES NOT WORK ON ASSEMBLIES)
    # -t <string> (input data type; HARDCODED TO '4' SINCE ONLY ASSEMBLIES ARE USED HERE)
    # -i <file> (/path/to/input/file)
    # -p <int> (number of threads for allele mode, if p > 4, only 4 threads will be used for assembly since the amount of extracted reads is small) 
    # -d <string> (output directory name)
    # -n <string> (sets the sample name in the report output)

    SeqSero2S.py \
      -m k \
      -p ~{cpu} \
      -t 4 \
      -n ~{samplename} \
      -d ~{samplename}_seqsero2s_output_dir \
      -i ~{assembly_fasta}

    # Run a python block to parse output file for terra data tables
    # contamination is not predicted when assembly data is used
    python3 <<CODE
    import csv
    with open("./~{samplename}_seqsero2s_output_dir/SeqSero_result.tsv",'r') as tsv_file:
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
    File seqsero2s_report = "./~{samplename}_seqsero2s_output_dir/SeqSero_result.tsv"
    String seqsero2s_version = read_string("VERSION")
    String seqsero2s_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero2s_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String seqsero2s_note = read_string("NOTE")
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