version 1.0

task extract_flu_segments {
  input {
    File assembly_fasta
    String flu_type # options: "Type_A" "Type_B"
    String flu_subtype # options: "Victoria" "Yamagata" "H3N2" "H1N1" "H5N1"
    Int cpu = 1
    Int memory = 2
    Int disk_size = 100
    String docker = "ubuntu"
  }
  String assembly_name = basename(basename(basename(assembly_fasta, ".fasta"), ".fa"), ".fna")
  command <<<
    set -euo pipefail

    # create directory to store the fastas of each individual segment
    mkdir -p temp_segments/

    # split the assembly fasta into segments based on headers
    csplit -s -z -f temp_segments/seq_ -b "%d" ~{assembly_fasta} '/^>/' '{*}'

    # get the lengths of each segment
    for file in temp_segments/seq_*; do
      if [ -s "$file" ]; then
        # extract the sequence (skip header)
        seq_length=$(grep -v "^>" "$file" | tr -d '\n' | wc -c)
        echo -e "$file\t$seq_length" >> temp_segments/lengths.tsv
      fi
    done

    # make sure there are 8 segments, if not, exit with an error
    segment_count=$(wc -l < temp_segments/lengths.tsv)
    if [ "$segment_count" -ne 8 ]; then
      echo "ERROR: Expected 8 segments, but found $segment_count. Please check the input assembly."
      exit 1
    fi

    # sort by length in descending order (largest to smallest)
    sort -k2,2nr temp_segments/lengths.tsv > temp_segments/sorted_lengths.tsv

    # make array of segment names based on flu type and order from largest to smallest
    if [[ ~{flu_type} == "Type_A" ]]; then
      segment_array=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")
    elif [[ ~{flu_type} == "Type_B" ]]; then
      segment_array=("PB1" "PB2" "PA" "HA" "NP" "NA" "MP" "NS")
    else
      echo "ERROR: Invalid flu type specified. Please use 'Type_A' or 'Type_B'."
      exit 1
    fi

    # copy each segment to a new file with the corresponding name
    i=0
    while IFS=$'\t' read -r file seq_length; do
      echo "Processing file: $file as segment ${segment_array[i]} with length $seq_length"
      segment_name="~{assembly_name}_~{flu_subtype}_${segment_array[i]}.fasta"
      cp "$file" "${segment_name}"
      i=$((i+1))
    done < temp_segments/sorted_lengths.tsv

    # remove all newlines from the input multi FASTA file to create a single line FASTA file and remove all headers and create a new headerline
    grep -v "^>" ~{assembly_fasta} | tr -d '\n' | sed '1i >~{assembly_name}_concatenated' > ~{assembly_name}_concatenated.fasta
    echo "" >> ~{assembly_name}_concatenated.fasta
  >>>
  output {
    File concatenated_fasta = "~{assembly_name}_concatenated.fasta"
    File seg_ha_assembly = "${assembly_name}_${flu_subtype}_HA.fasta"
    File seg_na_assembly = "${assembly_name}_${flu_subtype}_NA.fasta"
    File seg_pa_assembly = "${assembly_name}_${flu_subtype}_PA.fasta"
    File seg_pb1_assembly = "${assembly_name}_${flu_subtype}_PB1.fasta"
    File seg_pb2_assembly = "${assembly_name}_${flu_subtype}_PB2.fasta"
    File seg_mp_assembly = "${assembly_name}_${flu_subtype}_MP.fasta"
    File seg_np_assembly = "${assembly_name}_${flu_subtype}_NP.fasta"
    File seg_ns_assembly = "${assembly_name}_${flu_subtype}_NS.fasta"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 0
  }
}