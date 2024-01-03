version 1.0

task cat_files {
  input {
    Array[File] files_to_cat
    Array[String]? samplenames
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Boolean skip_extra_headers = false
    Boolean concatenate_variants = false
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<

    file_array=(~{sep=' ' files_to_cat})
    file_array_len=$(echo "${#file_array[@]}")
    samplename_array=(~{sep=' ' samplenames})
    samplename_array_len=$(echo "${#samplename_array[@]}")
    
    touch ~{concatenated_file_name}

    # Ensure file, and samplename arrays are of equal length
    if [ "$file_array_len" -ne "$samplename_array_len" ]; then
      echo "File array (length: $file_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
      exit 1
    fi

    if ! ~{concatenate_variants} ; then
      # cat files one by one and store them in the concatenated_files file, samplename will not be added as a column
      for index in ${!file_array[@]}; do
        file=${file_array[$index]}
        if ! ~{skip_extra_headers} ; then # act as if the first line of all files is not a header
          cat ${file} >> ~{concatenated_file_name}
        else # you want to skip the first line of all files except the first one
          if [ $index == 0 ]; then # if its the first file, cat the entire thing
            cat ${file} >> ~{concatenated_file_name}
          else # otherwise, skip the first line
            tail -n +2 ${file} >> ~{concatenated_file_name}
          fi
        fi
      done
    else
      # cat files one by one and store them in the concatenated_files file, but with an additional column indicating samplename
      for index in ${!file_array[@]}; do
        file=${file_array[$index]}
        samplename=${samplename_array[$index]}
        # create a new column with "samplename" as the column name and the samplename as the column content, combine with rest of file
        if [ "$index" -eq "0" ]; then
          # if first cloumn, add header
          awk -v var=$samplename 'BEGIN{ FS = OFS = "," } { print (NR==1? "samplename" : var), $0 }' $file > file.tmp
          cat file.tmp >> ~{concatenated_file_name}_concatenated_snps.csv
        else
          tail -n +2 $file | awk -v var=$samplename 'BEGIN{ FS = OFS = "," } { print var, $0 }' > file.tmp
          cat file.tmp >> ~{concatenated_file_name}_concatenated_snps.csv  
        fi
      done
    fi
  >>>
  output {
    File concatenated_files = "~{concatenated_file_name}"
  }
  runtime {
    docker: "~{docker_image}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task cat_files_fasta {
  input {
    Array[File] files_to_cat
    Array[String] headers
    String concatenated_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_cat})
    headers_array=(~{sep=' ' headers})
    touch ~{concatenated_file_name}

    # cat files one by one and store them in the concatenated_files file
    for index in ${!file_array[@]}; do
      file=${file_array[$index]}
      header=${headers_array[$index]}
      # replace the original header with the new header using sed before concatenating
      awk 1 ${file} | sed "s/^>.*/>${header}/" >> ~{concatenated_file_name}
    done
  >>>
  output {
    File concatenated_files = "~{concatenated_file_name}"
  }
  runtime {
    docker: "~{docker_image}"
    memory: "8 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}

task zip_files {
  input {
    Array[File] files_to_zip
    String zipped_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_zip})
    mkdir ~{zipped_file_name}

    # move files oto a single directory before zipping
    for index in ${!file_array[@]}; do
      file=${file_array[$index]}
      mv ${file} ~{zipped_file_name}
    done
    
    zip -r ~{zipped_file_name}.zip ~{zipped_file_name}
   >>>
  output {
    File zipped_files = "~{zipped_file_name}.zip"
  }
  runtime {
      docker: "~{docker_image}"
      memory: "8 GB"
      cpu: 2
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}

task transfer_files {
  input {
    Array[String] files_to_transfer
    String target_bucket
    Int cpus = 4
    Int mem_size_gb = 8
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_path_array="~{sep=' ' files_to_transfer}"

    gsutil -m cp -n ${file_path_array[@]} ~{target_bucket}
    
    echo "transferred_files" > transferred_files.tsv
    gsutil -m ls ~{target_bucket} >> transferred_files.tsv        
   >>>
  output {
    File transferred_files = "transferred_files.tsv"
  }
  runtime {
      docker: "~{docker_image}"
      memory: "~{mem_size_gb} GB"
      cpu: cpus
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}