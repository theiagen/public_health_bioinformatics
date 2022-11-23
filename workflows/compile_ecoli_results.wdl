version 1.0

workflow compile_results {

  input {
    Array[String]    SRR_array
    Array[File]      serotypefinder_array
    Array[File]      abricate_array
    Array[File]      abricate_virfinder_array
    Array[File]      amrfinder_array
  }
  call compile_abricate {
    input:
      array_srr=SRR_array,
      array_abr=abricate_array
  }

  call compile_abricate as compile_abricate_virfinder {
    input:
      array_srr=SRR_array,
      array_abr=abricate_virfinder_array
  }

  call compile_amrfinder {
    input:
      array_srr=SRR_array,
      array_afp=amrfinder_array
  }

  call compile_serotypefinder {
    input:
      array_srr=SRR_array,
      array_stf=serotypefinder_array
  }

  output {
    File      compiled_serotypefinder_results=compile_serotypefinder.compiled_results
    File      compiled_abricate_results=compile_abricate.compiled_results
    File      compiled_abricate_virfinder_results=compile_abricate_virfinder.compiled_results
    File      compiled_amrfinderplus_results=compile_amrfinder.compiled_results
  }
}


task compile_abricate {
  input {
    Array[String]     array_srr
    Array[File]       array_abr
  }

  command <<<
    touch results.txt

    srr_array=(~{sep=' ' array_srr})
    abr_array=(~{sep=' ' array_abr})
    echo "I am here"

    for index in ${!srr_array[@]}; do
      SRR=${srr_array[$index]}
      file=${abr_array[$index]}
      echo "$index"
      echo "$SRR"
      echo "$file"

      while IFS= read -r result
      do
      printf "%s %s\n" "$SRR $result" >> results.txt
      done < <(grep -E 'fasta' "$file")

    done
  >>>

  output {
    File      compiled_results="results.txt"
  }

  runtime {
    docker:       "quay.io/staphb/abricate:1.0.0"
    memory:       "4 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}

task compile_amrfinder {
  input {
    Array[String]     array_srr
    Array[File]       array_afp
  }

  command <<<
    touch results.txt

    srr_array=(~{sep=' ' array_srr})
    afp_array=(~{sep=' ' array_afp})
    echo "I am here"

    for index in ${!srr_array[@]}; do
      SRR=${srr_array[$index]}
      file=${afp_array[$index]}
      echo "$index"
      echo "$SRR"
      echo "$file"

      while IFS= read -r result
      do
      printf "%s %s\n" "$SRR $result" >> results.txt
      done < <(grep -E 'contig' "$file")

    done
  >>>

  output {
    File      compiled_results="results.txt"
  }

  runtime {
    docker:       "quay.io/staphb/ncbi-amrfinderplus:3.8.28"
    memory:       "4 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}


task compile_serotypefinder {
  input {
    Array[String]     array_srr
    Array[File]       array_stf
  }

  command <<<
    touch results.txt

    srr_array=(~{sep=' ' array_srr})
    stf_array=(~{sep=' ' array_stf})
    echo "I am here"

    for index in ${!srr_array[@]}; do
      SRR=${srr_array[$index]}
      file=${stf_array[$index]}
      echo "$index"
      echo "$SRR"
      echo "$file"

      while IFS= read -r result
      do
      printf "%s %s\n" "$SRR $result" >> results.txt
      done < <(grep -E 'fliC|wzy|wzx' "$file")

    done
  >>>

  output {
    File      compiled_results="results.txt"
  }

  runtime {
    docker:       "quay.io/staphb/serotypefinder:1.1"
    memory:       "4 GB"
    cpu:          1
    disks:        "local-disk 100 SSD"
    preemptible:  0
  }
}
