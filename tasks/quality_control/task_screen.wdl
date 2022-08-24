version 1.0

task check_reads {
  input {
    File read1
    File read2
    Int min_reads
    Int min_basepairs
    Int min_genome_size
    Int max_genome_size
    Int min_coverage
    Int min_proportion
    Boolean skip_screen
  }
  command <<<
    flag="PASS"
    estimated_genome_size=0
    if [[ "~{skip_screen}" = "false" ]] ; then
      
      # set cat command based on compression
      if [[ "~{read1}" == *".gz" ]] ; then
        cat_reads="zcat"
      else
        cat_reads="cat"
      fi

      # check one: number of reads
      read1_num=`eval "$cat_reads ~{read1}" | awk '{s++}END{print s/4}'`
      read2_num=`eval "$cat_reads ~{read2}" | awk '{s++}END{print s/4}'`
      # awk '{s++}END{print s/4' counts the number of lines and divides them by 4
      # key assumption: in fastq there will be four lines per read
      # sometimes fastqs do not have 4 lines per read, so this might fail one day

      if [ "${read1_num}" -le "~{min_reads}" ] || [ "${read2_num}" -le "~{min_reads}" ]; then
        flag="FAIL; the number of reads is below the minimum of ~{min_reads}"
      else
        flag="PASS"
      fi

      # checks two and three: number of basepairs and proportion of sequence
      if [ "${flag}" = "PASS" ]; then
        # count number of basepairs
        # this only works if the fastq has 4 lines per read, so this might fail one day
        read1_bp=`eval "${cat_reads} ~{read1}" | paste - - - - | cut -f2 | tr -d '\n' | wc -c`
        read2_bp=`eval "${cat_reads} ~{read2}" | paste - - - - | cut -f2 | tr -d '\n' | wc -c`
        # paste - - - - (print 4 consecutive lines in one row, tab delimited)
        # cut -f2 print only the second column (the second line of the fastq 4-line)
        # tr -d '\n' removes line endings
        # wc -c counts characters

        # set proportion variables for easy comparison
        # removing the , 2) to make these integers instead of floats
        percent_read1=$(python3 -c "print(round(($read1_bp / $read2_bp)*100))")
        percent_read2=$(python3 -c "print(round(($read2_bp / $read1_bp)*100))")

        if [ "$percent_read1" -lt "~{min_proportion}" ] ; then
          flag="FAIL; more than 50 percent of the total sequence is found in R2 (BP: $read2_bp; PERCENT: $percent_read2) compared to R1 (BP: $read1_bp; PERCENT: $percent_read1)"
        elif [ "$percent_read2" -lt "~{min_proportion}" ] ; then
          flag="FAIL; more than 50 percent of the total sequence is found in R1 (BP: $read1_bp; PERCENT: $percent_read1) compared to R2 (BP: $read2_bp; PERCENT: $percent_read2)"
        else
          flag="PASS"
        fi

        if [ "$flag" = "PASS" ] ; then
          if [ "${read1_bp}" -le "~{min_basepairs}" ] || [ "${read2_bp}" -le "~{min_basepairs}" ] ; then
            flag="FAIL; the number of basepairs is below the minimum of ~{min_basepairs}"
          else
            flag="PASS"
          fi
        fi    
      fi

      #checks four and five: estimated genome size and coverage
      if [ "${flag}" = "PASS" ]; then
        # determine genome size
              
        # First Pass; assuming average depth
        mash sketch -o test -k 31 -m 3 -r ~{read1} ~{read2} > mash-output.txt 2>&1
        grep "Estimated genome size:" mash-output.txt | \
          awk '{if($4){printf("%d", $4)}} END {if (!NR) print "0"}' > genome_size_output
        grep "Estimated coverage:" mash-output.txt | \
          awk '{if($3){printf("%d", $3)}} END {if (!NR) print "0"}' > coverage_output
        rm -rf test.msh
        rm -rf mash-output.txt
        estimated_genome_size=`head -n1 genome_size_output`
        estimated_coverage=`head -n1 coverage_output`

        # Check if second pass is needed
        if [ ${estimated_genome_size} -gt "~{max_genome_size}" ] || [ ${estimated_genome_size} -lt "~{min_genome_size}" ] ; then
          # Probably high coverage, try increasing number of kmer copies to 10
          M="-m 10"
          if [ ${estimated_genome_size} -lt "~{min_genome_size}" ]; then
            # Probably low coverage, try decreasing the number of kmer copies to 1
            M="-m 1"
          fi
          mash sketch -o test -k 31 ${M} -r ~{read1} ~{read2} > mash-output.txt 2>&1
          grep "Estimated genome size:" mash-output.txt | \
            awk '{if($4){printf("%d", $4)}} END {if (!NR) print "0"}' > genome_size_output
          grep "Estimated coverage:" mash-output.txt | \
            awk '{if($3){printf("%d", $3)}} END {if (!NR) print "0"}' > coverage_output
          rm -rf test.msh
          rm -rf mash-output.txt
        fi
        
        estimated_genome_size=`head -n1 genome_size_output`
        estimated_coverage=`head -n1 coverage_output`

        if [ "${estimated_genome_size}" -ge "~{max_genome_size}" ] ; then
          flag="FAIL; the genome size is estimated to be larger than the maximum of ~{max_genome_size} bps"
        elif [ "${estimated_genome_size}" -le "~{min_genome_size}" ] ; then
          flag="FAIL; the genome size is estimated to be smaller than the minimum of ~{min_genome_size} bps"
        else
          flag="PASS"   
          if [ "${estimated_coverage}" -le "~{min_coverage}" ] ; then
            flag="FAIL; the estimated coverage is less than the minimum of ~{min_coverage}x"
          else
            flag="PASS"
            echo $estimated_genome_size | tee EST_GENOME_LENGTH
          fi 
        fi
      fi 
    fi 
    
    echo $flag | tee FLAG
    echo $estimated_genome_size | tee EST_GENOME_LENGTH
  >>>
  output {
    String read_screen = read_string("FLAG")
    Int est_genome_length = read_int("EST_GENOME_LENGTH")
  }
  runtime {
    docker: "quay.io/bactopia/gather_samples:2.0.2"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}


task check_reads_se {
  input {
    File read1
    Int min_reads
    Int min_basepairs
    Int min_genome_size
    Int max_genome_size 
    Int min_coverage
    Boolean skip_screen 
  }
  command <<<
    flag="PASS"
    estimated_genome_size=0
    if [[ "~{skip_screen}" = "false" ]] ; then
      
      # set cat command based on compression
      if [[ "~{read1}" == *".gz" ]] ; then
        cat_reads="zcat"
      else
        cat_reads="cat"
      fi

      # check one: number of reads
      read1_num=`eval "$cat_reads ~{read1}" | awk '{s++}END{print s/4}'`
      # awk '{s++}END{print s/4' counts the number of lines and divides them by 4
      # key assumption: in fastq there will be four lines per read
      # sometimes fastqs do not have 4 lines per read, so this might fail one day

      if [ "${read1_num}" -le "~{min_reads}" ] ; then
        flag="FAIL; the number of reads is below the minimum of ~{min_reads}"
      else
        flag="PASS"
      fi

      # checks two and three: number of basepairs and proportion of sequence
      if [ "${flag}" = "PASS" ]; then
        # count number of basepairs
        # this only works if the fastq has 4 lines per read, so this might fail one day
        read1_bp=`eval "${cat_reads} ~{read1}" | paste - - - - | cut -f2 | tr -d '\n' | wc -c`
        # paste - - - - (print 4 consecutive lines in one row, tab delimited)
        # cut -f2 print only the second column (the second line of the fastq 4-line)
        # tr -d '\n' removes line endings
        # wc -c counts characters

        if [ "$flag" = "PASS" ] ; then
          if [ "${read1_bp}" -le "~{min_basepairs}" ] ; then
            flag="FAIL; the number of basepairs is below the minimum of ~{min_basepairs}"
          else
            flag="PASS"
          fi
        fi    
      fi

      #checks four and five: estimated genome size and coverage
      if [ "${flag}" = "PASS" ]; then
        # determine genome size
              
        # First Pass; assuming average depth
        mash sketch -o test -k 31 -m 3 -r ~{read1} > mash-output.txt 2>&1
        grep "Estimated genome size:" mash-output.txt | \
          awk '{if($4){printf("%d", $4)}} END {if (!NR) print "0"}' > genome_size_output
        grep "Estimated coverage:" mash-output.txt | \
          awk '{if($3){printf("%d", $3)}} END {if (!NR) print "0"}' > coverage_output
        rm -rf test.msh
        rm -rf mash-output.txt
        estimated_genome_size=`head -n1 genome_size_output`
        estimated_coverage=`head -n1 coverage_output`

        # Check if second pass is needed
        if [ ${estimated_genome_size} -gt "~{max_genome_size}" ] || [ ${estimated_genome_size} -lt "~{min_genome_size}" ] ; then
          # Probably high coverage, try increasing number of kmer copies to 10
          M="-m 10"
          if [ ${estimated_genome_size} -lt "~{min_genome_size}" ]; then
            # Probably low coverage, try decreasing the number of kmer copies to 1
            M="-m 1"
          fi
          mash sketch -o test -k 31 ${M} -r ~{read1} > mash-output.txt 2>&1
          grep "Estimated genome size:" mash-output.txt | \
            awk '{if($4){printf("%d", $4)}} END {if (!NR) print "0"}' > genome_size_output
          grep "Estimated coverage:" mash-output.txt | \
            awk '{if($3){printf("%d", $3)}} END {if (!NR) print "0"}' > coverage_output
          rm -rf test.msh
          rm -rf mash-output.txt
        fi
        
        estimated_genome_size=`head -n1 genome_size_output`
        estimated_coverage=`head -n1 coverage_output`

        if [ "${estimated_genome_size}" -ge "~{max_genome_size}" ] ; then
          flag="FAIL; the genome size is estimated to be larger than the maximum of ~{max_genome_size} bps"
        elif [ "${estimated_genome_size}" -le "~{min_genome_size}" ] ; then
          flag="FAIL; the genome size is estimated to be smaller than the minimum of ~{min_genome_size} bps"
        else
          flag="PASS"   
          if [ "${estimated_coverage}" -le "~{min_coverage}" ] ; then
            flag="FAIL; the estimated coverage is less than the minimum of ~{min_coverage}x"
          else
            flag="PASS"
            echo $estimated_genome_size | tee EST_GENOME_LENGTH
          fi 
        fi
      fi 
    fi 
    
    echo $flag | tee FLAG
    echo $estimated_genome_size | tee EST_GENOME_LENGTH
  >>>
  output {
    String read_screen = read_string("FLAG")
    Int est_genome_length = read_int("EST_GENOME_LENGTH")
  }
  runtime {
    docker: "quay.io/bactopia/gather_samples:2.0.2"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
    maxRetries: 3
  }
}