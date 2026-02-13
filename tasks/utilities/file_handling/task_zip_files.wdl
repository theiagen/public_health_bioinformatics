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
    
      echo "DEBUG: Pulling $file"
      if [ -f "$file" ]; then
        echo "DEBUG: $file exists"
        filename=$(basename "$file") # Extract the filename (e.g., test.tsv)
        dest="~{zipped_file_name}/$filename"

        # Counter is always set to 1 so that if there are 
        # other duplicated filenames they will be counted as well.
        counter=1

        echo "DEBUG: Checking for $file in $dest"
        # Check for duplicate files in the destination
        while [ -e "$dest" ]; do 
          echo "DEBUG: Duplicate filename found, adding a file index for differentiation."
          dest="~{zipped_file_name}/${filename%.*}_${counter}.${filename##*.}"
          echo "DEBUG: New filename ${filename%.*}_${counter}.${filename##*.}"
          ((counter++))
        done

        # Move the file to the destination with the new name
        # If loop is not entered, filename will remain unchanged. 
        mv "$file" "$dest"

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