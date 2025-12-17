version 1.0

task genoflu {
  input {
    File assembly_fasta
    String samplename
    Float min_percent_identity = 98 # genoflu default is 98
    # excel file to cross-reference BLAST findings; probably useful if novel
    #  genotypes are not in the default file used by genoflu.py
    File? cross_reference
    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/genoflu:1.06"
    Int memory = 2
  }
  command <<<
    set -euo pipefail

    cp ~{assembly_fasta} .

    echo "DEBUG: capturing genoflu version..."
    genoflu.py -v | sed -e 's/genoflu.py:\ version\ //' | tee VERSION

    echo "DEBUG: running genoflu.py..."
    genoflu.py \
      --fasta ~{assembly_fasta} \
      --sample_name ~{samplename} \
      ~{"--pident_threshold " + min_percent_identity} \
      ~{"--cross_reference" + cross_reference} > genoflu.output.txt

    GENOTYPE=$(grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2 | cut -d " " -f2 | cut -d ":" -f1)
    ALL_SEGMENTS=$(grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2 | cut -d " " -f3-)

    # If genotype unable to be assigned ("Not"), then parse out the expected text
    if [[ "$GENOTYPE" == "Not" ]]; then
      echo "DEBUG: parsing out genotype and all segments..."
      grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2- | cut -d " " -f2- | cut -d ":" -f1 | tee GENOTYPE
      grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2- | cut -d " " -f4- | cut -d ":" -f1 | tee ALL_SEGMENTS
    else
      echo "DEBUG: saving genotype and all segments..."
      echo "$GENOTYPE" | tee GENOTYPE
      echo "$ALL_SEGMENTS" | tee ALL_SEGMENTS
    fi
    
    mv -v ~{samplename}_*_stats.tsv ~{samplename}_stats.tsv
  >>>
  output {
    String genoflu_version = read_string("VERSION")
    String genoflu_genotype = read_string("GENOTYPE")
    String genoflu_all_segments = read_string("ALL_SEGMENTS")
    File genoflu_output_tsv = "~{samplename}_stats.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}