version 1.0

task zip_files {
  input {
    Array[File] files_to_zip
    String zipped_file_name
    String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
    Int memory = 8
    Int cpu = 2
    Int disk_size = 100
  }
  meta {
    # added so that call caching is always turned off
    volatile: true
  }
  command <<<
    file_array=(~{sep=' ' files_to_zip})
    mkdir ~{zipped_file_name}

    # move files into a single directory before zipping
    for file in "${file_array[@]}"; do
      echo "debug file found $file"
      if [ -f "$file" ]; then
        echo "debug file found $file"
        filename=$(basename "$file") # Extract the filename (e.g., test.tsv)
        ext="${filename##*.}"        # Extract the file extension (e.g., tsv)
        name="${filename%.*}"        # Extract the file name without extension (e.g., test)

        echo "DEBUG filename $filename"
        echo "DEBUG ext $ext"
        echo "DEBUG name $name"
        
        # Counter is always set to 1 so that if there are 
        # other duplicated filenames they will be counted as well.
        counter=1
        new_name="$name.$ext"

        # Check for duplicate files in the destination
        while [ -e "~{zipped_file_name}/$new_name" ]; do 
          echo "DEBUG file exists"
          new_name="${name}_${counter}.$ext"
          ((counter++))
        done

        # Move the file to the destination with the new name
        mv "$file" "~{zipped_file_name}/$new_name"
      else
        echo "File not found: $file"
      fi
    done
    
    zip -r ~{zipped_file_name}.zip ~{zipped_file_name}
   >>>
  output {
    File zipped_files = "~{zipped_file_name}.zip"
  }
  runtime {
    docker: "~{docker_image}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}