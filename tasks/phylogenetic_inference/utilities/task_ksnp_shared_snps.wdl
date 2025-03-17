version 1.0

task ksnp_shared_snps {
  input {
    File ksnp_vcf_ref_genome
    Array[String] samplename
    String cluster_name
    Int disk_size = 100
  }
  command <<<

    python3 <<CODE
    import pandas as pd

    file_path="~{ksnp_vcf_ref_genome}"
    samplename_array = "~{sep=',' samplename}"

    # split the samplename_array into a list 
    # this will be used to isolate only columns in the df with sample data
    sample_columns = samplename_array.split(',')

    # we want to skip the commented lines and use the line starting with #CHROM as the header
    # grab the line that starts with #CHROM 
    header = None
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                break

    # Read the VCF file skipping comment lines and specifying the header
    data = pd.read_csv(file_path, sep='\t', comment='#', header=None, names=header)

    # Function to count occurrences of '.'
    # a '.' indicates that the same does not have the k-mer, so it is not a core k-mer
    def count_dot(row):
      return sum(1 for value in row if value == '.')

    # Apply the function to each row to count occurrences of '.'
    data['num_dots'] = data[sample_columns].apply(count_dot, axis=1)

    # Drop rows where the number of '.' is greater than 0
    data_core = data[data['num_dots'] == 0].drop(columns='num_dots')

    # Convert the values in selected columns to numeric 
    data_core[sample_columns] = data_core[sample_columns].apply(pd.to_numeric, errors='coerce')

    # Sum values by row after the 'FORMAT' column
    data_core['Row_Sum'] = data_core[sample_columns].sum(axis=1)

    # Sort DataFrame by the 'Row_Sum' column in descending order
    sorted_data_core = data_core.sort_values(by='Row_Sum', ascending=False).drop(columns='Row_Sum')

    # Write the sorted DataFrame to a CSV file
    sorted_data_core.to_csv('~{cluster_name}_core_snp_table.csv', index=False)

    CODE
  >>>
  output{
    File core_snp_table = "~{cluster_name}_core_snp_table.csv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/staphb/mykrobe:0.12.1" # used because it contains both biopython and pandas
    memory: "2 GB"
    cpu: 2
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}