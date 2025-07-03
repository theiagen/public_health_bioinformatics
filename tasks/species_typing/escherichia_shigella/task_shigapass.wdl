version 1.0

task shigapass_many {
  meta {
    description: "In-silico prediction of Shigella serotypes & EIEC differentiation. This task is designed to run ShigaPass on multiple assemblies at once, and requires an array of assemblies and sample names."
  }
  input {
    Array[File] assembly_fasta
    Array[String] samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/shigapass:1.5.0"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 8
  }
  command <<<
    # prevent failure; especially important since ShigaPass is a shell script
    set -euo pipefail
    # capture version
    ShigaPass.sh -v | sed 's/ShigaPass version //g' | tee VERSION

    assembly_array=(~{sep=' ' assembly_fasta})
    assembly_array_len=$(echo "${#assembly_array[@]}")
    samplename_array=(~{sep=' ' samplename})
    samplename_array_len=$(echo "${#samplename_array[@]}")

    # Ensure assembly, and samplename arrays are of equal length
    if [ "$assembly_array_len" -ne "$samplename_array_len" ]; then
      echo "Assembly array (length: $assembly_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
      exit 1
    fi

    # create a txt file of the assembly paths
    echo "Creating ASSEMBLY_PATH file..."
    for index in "${!assembly_array[@]}"; do
      assembly=${assembly_array[$index]}
      realpath "${assembly}" | tee -a ASSEMBLY_PATHS
    done
    echo "Finished creating ASSEMBLY_PATH file. Contents:"
    cat ASSEMBLY_PATHS
    echo 

    # run ShigaPass; FYI database path input is mandatory and is specific to this docker image
    echo "running ShigaPass..."
    ShigaPass.sh \
      -l ASSEMBLY_PATHS \
      -p /ShigaPass-1.5.0/SCRIPT/ShigaPass_DataBases/ \
      -t ~{cpu} \
      -k \
      -o shigapass/

    echo "Finished running ShigaPass."

    # convert summary file from semi-colon separate values into TSV format
    echo "Converting summary file to TSV format..."
    sed 's/;/\t/g' shigapass/ShigaPass_summary.csv > shigapass/ShigaPass_summary.tsv

    # renaming file to end in .txt since it is semicolon separated values sheet. Easier to open in Excel if ends in .txt
    mv -v shigapass/ShigaPass_summary.csv shigapass/ShigaPass_summary.txt

    # if the ShigaPass_Flex_summary.csv file exists, convert to tabular format and rename the original file to end in .txt
    if [ -f shigapass/ShigaPass_Flex_summary.csv ]; then
      # convert to TSV
      sed 's/;/\t/g' shigapass/ShigaPass_Flex_summary.csv > shigapass/ShigaPass_Flex_summary.tsv

      echo "Renaming ShigaPass_Flex_summary.csv to ShigaPass_Flex_summary.txt"
      mv -v shigapass/ShigaPass_Flex_summary.csv shigapass/ShigaPass_Flex_summary.txt
    else
      echo "ShigaPass_Flex_summary.csv does not exist, skipping renaming."
    fi

    # rename intermediate files to include subdir name to prevent duplicate file names in glob output
    echo "Renaming intermediate files to include subdir name..."
    for file in shigapass/*/*.txt; do
      # get the base name of the file (e.g., ShigaPass_Flex_summary.txt)
      base_name=$(basename "$file") 
      # get the directory name (e.g., shigapass/ShigaPass_Flex)
      dir_name=$(dirname "$file")
      # create the new file name with the subdir name included
      new_file_name="${dir_name}/$(basename "$dir_name")_${base_name%.txt}.txt"
      # rename the file
      mv -v "$file" "$new_file_name"
    done 
  >>>
  output {
    File shigapass_summary = "shigapass/ShigaPass_summary.txt"
    File shigapass_summary_tsv = "shigapass/ShigaPass_summary.tsv"
    File? shigapass_flexneri_summary = "shigapass/ShigaPass_Flex_summary.txt"
    File? shigapass_flexneri_summary_tsv = "shigapass/ShigaPass_Flex_summary.tsv"
    # keeping the intermediate files in the output for now, but may remove later. Not planning to output at workflow level.
    # NOTE: the intermediate files are renamed to include the subdir name to prevent duplicate file names
    Array[File] intermediate_files = glob("shigapass/*/*.txt")
    String shigapass_version = read_string("VERSION")
    String shigapass_docker = docker
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    # keep at 0 since ShigaPass can take a long time on large batches (>100 samples)
    preemptible: 0
  }
}
