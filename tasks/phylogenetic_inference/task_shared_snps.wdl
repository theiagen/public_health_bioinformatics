version 1.0

task shared_snps {
  input {
    Array[File] snippy_variants_results
    Array[String] samplenames
    String concatenated_file_name

    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    Int disk_size = 100
    Int cpu = 1
    Int memory = 8
  }
  command <<<
    
    file_array=(~{sep=" " snippy_variants_results})
    file_array_len=$(echo "${#file_array[@]}")
    samplename_array=(~{sep=' ' samplenames})
    samplename_array_len=$(echo "${#samplename_array[@]}")
    
    touch ~{concatenated_file_name}_concatenated_snps.csv

    # Ensure file, and samplename arrays are of equal length
    if [ "$file_array_len" -ne "$samplename_array_len" ]; then
      echo "File array (length: $file_array_len) and samplename array (length: $samplename_array_len) are of unequal length." >&2
      exit 1
    fi

    # cat files one by one and store them in the concatenated_files file

    for index in ${!file_array[@]}; do
      if [ "$index" -eq "0" ]; then
        file=${file_array[$index]}
        samplename=${samplename_array[$index]}
        # create a new column with "samplename" as the column name and the samplename as the column content, combine with rest of file
        awk -v var=$samplename 'BEGIN{ FS = OFS = "," } { print (NR==1? "samplename" : var), $0 }' $file > file.tmp
        cat file.tmp >> ~{concatenated_file_name}_concatenated_snps.csv
      else
        file=${file_array[$index]}
        samplename=${samplename_array[$index]}
        tail -n +2 $file | awk -v var=$samplename 'BEGIN{ FS = OFS = "," } { print var, $0 }' > file.tmp
        cat file.tmp >> ~{concatenated_file_name}_concatenated_snps.csv  
      fi
    done

    python3 <<CODE 
    import pandas as pd

    # read the concatenated file into terra
    tablename = "~{concatenated_file_name}_concatenated_snps.csv" 
    df = pd.read_csv(tablename, delimiter=',', header=0, index_col=False) 
    
    # fill empty columns with NA
    df.fillna("NA", inplace=True)

    # parse evidence column?

    # pivot the table so that the samples are columns and 
    # CHROM POS LOCUS_TAG PRODUCT (concatenated) are rows
    # in the process, fill the pivoted table with the ALT column 
    # and fill with 0 if there is no ALT
    pivoted_df = df.pivot_table(index=['CHROM','POS','TYPE','REF','ALT','FTYPE','STRAND','NT_POS','AA_POS','EFFECT','LOCUS_TAG','GENE','PRODUCT'], columns='samplename', values='EVIDENCE', aggfunc='first')
    
    # fill with 0 if there is no ALT
    pivoted_df = pivoted_df.fillna(0)

    # Replace non-zero values with 1
    pivoted_df_binary = pivoted_df.applymap(lambda x: 1 if x != 0 else 0)

    # write as csv
    pivoted_df_binary.to_csv('~{concatenated_file_name}_shared_snp_table.csv')

    CODE
  >>>
  output {
    File snippy_concatenated_snps = "~{concatenated_file_name}_concatenated_snps.csv"
    File snippy_shared_snp_table = "~{concatenated_file_name}_shared_snp_table.csv" 
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}