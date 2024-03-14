version 1.0

task shared_variants {
  input {
    File concatenated_variants
    String concatenated_file_name
    String docker = "us-docker.pkg.dev/general-theiagen/theiagen/terra-tools:2023-03-16"
    Int disk_size = 100
    Int cpu = 1
    Int memory = 8
  }
  command <<<
  
    python3 <<CODE 
    import pandas as pd

    # read the concatenated file into terra
    tablename = "~{concatenated_variants}" 
    df = pd.read_csv(tablename, delimiter=',', header=0, index_col=False) 
    
    # fill empty columns with NA (particularly relevant for mutations with no associated AA changes)
    df.fillna("NA", inplace=True)

    # pivot the table so that the samples are columns and 
    # CHROM,POS,TYPE,REF,ALT,FTYPE,STRAND,NT_POS,AA_POS,EFFECT,LOCUS_TAG,GENE,PRODUCT are rows
    # in the process, fill the pivoted table with the EVIDENCE column 
    pivoted_df = df.pivot_table(index=['CHROM','POS','TYPE','REF','ALT','FTYPE','STRAND','NT_POS','AA_POS','EFFECT','LOCUS_TAG','GENE','PRODUCT'], columns='samplename', values='EVIDENCE', aggfunc='first')
    
    # if a sample does not have a mutation that is present in another sample, 
    # the pivot table will have a NaN value
    # so replace NaN with 0
    pivoted_df = pivoted_df.fillna(0)

    # Replace non-zero values with 1
    pivoted_df_binary = pivoted_df.applymap(lambda x: 1 if x != 0 else 0)

    # write as csv
    pivoted_df_binary.to_csv('~{concatenated_file_name}_shared_snp_table.csv')

    CODE
  >>>
  output {
    File shared_variants_table = "~{concatenated_file_name}_shared_snp_table.csv" 
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