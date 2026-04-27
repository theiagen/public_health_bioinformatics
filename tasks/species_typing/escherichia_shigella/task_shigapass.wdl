version 1.0

task shigapass {
  input {
    File assembly
    String samplename
    
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/shigapass:1.5.0"
    Int disk_size = 50
    Int cpu = 2
    Int memory = 4
  }
  command <<<
    set -euo pipefail
    
    ShigaPass.sh -v | sed 's/ShigaPass version //g' | tee VERSION
    
    # shigapass requires a text file containing the filepath instead of the filepath itself 
    realpath ~{assembly} > assembly_path.txt

    ShigaPass.sh \
      -l assembly_path.txt \
      -p /ShigaPass-1.5.0/SCRIPT/ShigaPass_DataBases/ \
      -t ~{cpu} \
      -k \
      -o shigapass/
   
    # convert summary files from semi-colon separate values into TSV format
    sed 's/;/\t/g' shigapass/ShigaPass_summary.csv > shigapass/ShigaPass_summary.tsv

    if [ -f shigapass/ShigaPass_Flex_summary.csv ]; then
      sed 's/;/\t/g' shigapass/ShigaPass_Flex_summary.csv > shigapass/ShigaPass_Flex_summary.tsv
    else
      echo "ShigaPass_Flex_summary.csv does not exist, skipping renaming."
    fi

    cut -f6 -d';' shigapass/ShigaPass_summary.csv | tail -n1 > IPAH_PRESENCE_ABSENCE
    cut -f7 -d';' shigapass/ShigaPass_summary.csv | tail -n1 > PREDICTED_SEROTYPE
    cut -f9 -d';' shigapass/ShigaPass_summary.csv | tail -n1 > PREDICTED_SEROTYPE_FLEXNERI
  >>>
  output {
    String shigapass_predicted_serotype = read_string("PREDICTED_SEROTYPE")
    String shigapass_predicted_serotype_flexneri = read_string("PREDICTED_SEROTYPE_FLEXNERI")
    String shigapass_ipaH_presence_absence = read_string("IPAH_PRESENCE_ABSENCE")

    File shigapass_summary_tsv = "shigapass/ShigaPass_summary.tsv"
    File? shigapass_flexneri_summary_tsv = "shigapass/ShigaPass_Flex_summary.tsv"
    Array[File] shigapass_intermediate_files = glob("~{samplename}/*")

    String shigapass_version = read_string("VERSION")
    String shigapass_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}