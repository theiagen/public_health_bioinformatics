version 1.0

task consensus_qc {
  input {
    File assembly_fasta
    File? reference_genome
    Int? genome_length
    Int disk_size = 100
  }
  command <<<
    if [ ~{reference_genome} ] ; then
      GENOME_LEN=$(grep -v ">" ~{reference_genome} | tr --delete '\n' | wc -c)
    elif [ ~{genome_length} ] ; then
      GENOME_LEN=~{genome_length}
    else
      # set SC2 default
      GENOME_LEN=29903
    fi 

    # capture date and version
    date | tee DATE

    num_N=$( grep -v ">" ~{assembly_fasta} | grep -o 'N' | wc -l )
    if [ -z "$num_N" ] ; then num_N="0" ; fi
    echo $num_N | tee NUM_N

    num_ACTG=$( grep -v ">" ~{assembly_fasta} | grep -o -E "C|A|T|G" | wc -l )
    if [ -z "$num_ACTG" ] ; then num_ACTG="0" ; fi
    echo $num_ACTG | tee NUM_ACTG

    # calculate percent coverage (Wu Han-1 genome length: 29903bp)
    python3 -c "print ( round( ($num_ACTG / $GENOME_LEN ) * 100, 2 ) )" | tee PERCENT_REF_COVERAGE

    num_degenerate=$( grep -v ">" ~{assembly_fasta} | grep -o -E "B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z" | wc -l )
    if [ -z "$num_degenerate" ] ; then num_degenerate="0" ; fi
    echo $num_degenerate | tee NUM_DEGENERATE

    num_total=$( grep -v ">" ~{assembly_fasta} | grep -o -E '[A-Z]' | wc -l )
    if [ -z "$num_total" ] ; then num_total="0" ; fi
    echo $num_total | tee NUM_TOTAL
  >>>
  output {
    Int number_N = read_string("NUM_N")
    Int number_ATCG = read_string("NUM_ACTG")
    Int number_Degenerate = read_string("NUM_DEGENERATE")
    Int number_Total = read_string("NUM_TOTAL")
    Float percent_reference_coverage = read_string("PERCENT_REF_COVERAGE")
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    memory: "2 GB"
    cpu: 1
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" # TES
    preemptible: 0
  }
}