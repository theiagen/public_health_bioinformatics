version 1.0

task emmtypingtool {
  meta {
    description: "emm-typing of Streptococcus pyogenes raw reads"
  }
  input {
    File read1
    File read2
    String samplename
    String docker = "staphb/emmtypingtool:0.0.1"
    Int cpu = 2

    # Parameters
    #-input_directory INPUT_DIRECTORY, -i INPUT_DIRECTORY
    #--fastq_1 FASTQ_1, -1 FASTQ_1
    #--fastq_2 FASTQ_2, -2 FASTQ_2
    #--profile_file_directory PROFILE_FILE_DIRECTORY, -m PROFILE_FILE_DIRECTORY
    #--output_directory OUTPUT_DIRECTORY, -o OUTPUT_DIRECTORY
    #--bowtie BOWTIE, -b BOWTIE
    #--samtools SAMTOOLS, -sam SAMTOOLS
    #--log_directory LOG_DIRECTORY, -log LOG_DIRECTORY
  }
  command <<<
    #echo $(emmtyper --version 2>&1) | sed 's/^.*emmtyper v//' | tee VERSION
    emm_typing.py \
      -m /db \
      -1 ~{read1} \
      -2 ~{read2} \
      -o output_dir

    grep "version" output_dir/*.results.xml | sed -n 's/.*version="\([^"]*\)".*/\1/p' | tee VERSION
    grep "Final_EMM_type" output_dir/*.results.xml | sed -n 's/.*value="\([^"]*\)".*/\1/p' | tee EMM_type
    mv output_dir/*.results.xml ~{samplename}.xml
  >>>
  output {
    String emmtypingtool_results = read_string("EMM_type")
    File emmtypingtool_xml = "~{samplename}.xml"
    String emmtypingtool_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk 50 SSD"
    preemptible: 0
  }
}
