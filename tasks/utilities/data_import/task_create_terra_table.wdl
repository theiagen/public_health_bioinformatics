version 1.0

task create_terra_table {
  input {
    String new_table_name
    String data_location_path
    String? file_ending # comma-delimited list
    Boolean paired_end
    Boolean assembly_data

    String terra_project
    String terra_workspace

    Int disk_size = 25
    Int cpu = 1
    Int memory = 4
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-06-21"
  }
  meta {
    volatile: true
  }
  String new_table_name_updated = sub(new_table_name, " ", "_")
  command <<<
    echo "DEBUG: starting to creating a terra table"

    # writing these to stderr so if errors occur, the user can check what they provided.
    echo "DEBUG: new_table_name: ~{new_table_name_updated}" >&2
    echo "DEBUG: data_location_path: ~{data_location_path}" >&2
    echo "DEBUG: paired_end: ~{paired_end}" >&2
    echo "DEBUG: assembly_data: ~{assembly_data}" >&2 

    if ~{paired_end} && ~{assembly_data}; then
      echo "DEBUG: paired-end AND assembly data indicated; this is not supported. Please check your input parameters and then try again" >&2
      exit 1
    elif ~{assembly_data}; then
      echo "DEBUG: single-end (paired-end = false) AND assembly data indicated; only assembly data will be retrieved." >&2
    fi

    if ~{paired_end}; then
      echo "DEBUG: paired-end data indicated"
      #  this pattern matches any files with either _R1, _1, _R2, or _2 in the filename that ends in .fastq(.gz)
      PATTERN="_R*[1-2].*\b\.f(q|astq)(\.gz)?\b$"
      echo -e "entity:~{new_table_name_updated}_id\tread1\tread2\tupload_date\ttable_created_by" > terra_table_to_upload.tsv
    else
      echo "DEBUG: single-end data indicated"
      # this pattern matches any files that end in fastq(.gz)
      PATTERN="\b\.f(q|astq)(\.gz)?\b$"
      echo -e "entity:~{new_table_name_updated}_id\tread1\tupload_date\ttable_created_by" > terra_table_to_upload.tsv
    fi

    if ~{assembly_data}; then
      echo "DEBUG: assembly data indicated"
      # this pattern matches any files that end in .fasta(.gz), .fa(.gz) or .fna(.gz)
      PATTERN="\b\.f(na|a|as|fn|asta)(.gz)?\b$"
      echo -e "entity:~{new_table_name_updated}_id\tassembly_fasta\tupload_date\ttable_created_by" > terra_table_to_upload.tsv
    fi

    if [ -n "~{file_ending}" ]; then
      echo "DEBUG: a file ending pattern was provided; using that instead of the default pattern"
      # this could be a comma delimited list, so we will want to split that up
      IFS="," read -ra FILE_ENDINGS <<< "~{file_ending}"

      echo ${FILE_ENDINGS[@]} >&2
      
      PATTERN="("
      for ending in "${FILE_ENDINGS[@]}"; do
        if [ $PATTERN == "(" ]; then
          PATTERN="$PATTERN$ending"
        else
         PATTERN="$PATTERN|$ending"
        fi
      done

      PATTERN="$PATTERN)"
      echo "$PATTERN"
      unset IFS
    fi

    # get list of files
    gcloud storage ls ~{data_location_path} | grep -E "$PATTERN" > filelist-fullpath.txt

    if [ -s filelist-fullpath.txt ]; then
      echo "DEBUG: files were identified; now identifying sample names"
    else
      echo "DEBUG: no files were identified; please check your data location path and then try again" >&2
      exit 1
    fi

    # use basename to prevent grep searches matching with file paths that include a _(R)1 or _(R)2
    while read -r filepath; do
      basename "$filepath" >> filelist-filename.txt
    done <filelist-fullpath.txt

    UPLOAD_DATE=$(date -I)

    touch samplenames.txt
    while read -r filepath; do
      file=$(basename "$filepath")

      if [ -n "~{file_ending}" ]; then
        echo "DEBUG: file ending pattern was provided; using that instead of the default pattern"
        samplename=$file
        for ending in "${FILE_ENDINGS[@]}"; do
          samplename=${samplename%%$ending}
        done
      else
        echo "DEBUG: no file ending pattern was provided; using the default pattern"
        # samplename is everything before the first underscore and first decimal (name.banana.hello_yes_please.fastq.gz -> name)
        no_underscore_samplename=${file%%_*} 
        samplename=${no_underscore_samplename%%.*}
      fi


      if grep "$samplename" samplenames.txt; then
        echo "DEBUG: $samplename is in the terra table"
      else
        echo "DEBUG: $samplename has NOT been added to the terra table"
        echo "$samplename" >> samplenames.txt

        if ~{paired_end}; then
          READ1_PATTERN="_R*1.*\b\.f(q|astq)(\.gz)?\b$"
          READ2_PATTERN="_R*2.*\b\.f(q|astq)(\.gz)?\b$"

          # search for the appropriate file in the list of filenames that exclude the path (filelist-filename.txt) 
          #  and then search for that file in the full-path filelist (filelist-fullpath.txt)
          read1=$(grep $(grep -E "$READ1_PATTERN" filelist-filename.txt | grep "$samplename") filelist-fullpath.txt)
          read2=$(grep $(grep -E "$READ2_PATTERN" filelist-filename.txt | grep "$samplename") filelist-fullpath.txt)

          # occasionally the readnames will be prefixed with `filelist-fullpath.txt:` so we need to remove that
          read1=$(echo $read1 | sed -e 's/^filelist-fullpath.txt://')
          read2=$(echo $read2 | sed -e 's/^filelist-fullpath.txt://')
          
          echo -e "$samplename\t$read1\t$read2\t$UPLOAD_DATE\tCreate_Terra_Table_PHB" >> terra_table_to_upload.tsv 
        else
          echo -e "$samplename\t$filepath\t$UPLOAD_DATE\tCreate_Terra_Table_PHB" >> terra_table_to_upload.tsv
        fi

      fi
    done <filelist-fullpath.txt

    echo "DEBUG: terra table created, now beginning upload"
    python3 /scripts/import_large_tsv/import_large_tsv.py --project "~{terra_project}" --workspace "~{terra_workspace}" --tsv terra_table_to_upload.tsv
  >>>
  output {
    File terra_table_to_upload = "terra_table_to_upload.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 1
  }
}