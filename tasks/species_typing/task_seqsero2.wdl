version 1.0

task seqsero2 {
  # Inputs
  input {
    File read1
    File? read2
    String samplename
    String mode ="a"
    String seqsero2_docker_image = "quay.io/staphb/seqsero2:1.2.1"
    Boolean paired_end
  }

  command <<<
    # capture date and version
    # Print and save date
    date | tee DATE
    # Print and save version
    SeqSero2_package.py --version | tee VERSION
    # Run SeqSero2 on the input read data
    SeqSero2_package.py \
    -p 8 \
    ~{true='-t 2' false='-t 3' paired_end} \
    -m ~{mode} \
    -n ~{samplename} \
    -d ~{samplename}_seqseqro2_output_dir \
    -i ~{read1} ~{read2}
    # Run a python block to parse output file for terra data tables
    python3 <<CODE
    import csv
    with open("./~{samplename}_seqseqro2_output_dir/SeqSero_result.tsv",'r') as tsv_file:
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

    CODE
  >>>
  output {
    File seqsero2_report = "./~{samplename}_seqseqro2_output_dir/SeqSero_result.tsv"
    String seqsero2_version = read_string("VERSION")
    String seqsero2_predicted_antigenic_profile = read_string("PREDICTED_ANTIGENIC_PROFILE")
    String seqsero2_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String seqsero2_predicted_contamination = read_string("CONTAMINATION")
  }
  runtime {
    docker:       "~{seqsero2_docker_image}"
    memory:       "16 GB"
    cpu:          8
    disks:        "local-disk 100 SSD"
    preemptible:  0
    maxRetries:   3
  }
}