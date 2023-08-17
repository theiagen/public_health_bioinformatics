version 1.0

task tb_gene_coverage {
  input {
    File bamfile
    File bamindex
    String samplename
    Int min_depth = 10
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
  }
  command <<<

    # grab file from tbdb with reference positions for each resistance gene
    # source: https://github.com/jodyphelan/TBProfiler/blob/master/db/tbdb.bed
    wget https://raw.githubusercontent.com/jodyphelan/TBProfiler/master/db/tbdb.bed

    # parse file and output coverage with samtools depth
    python3 <<CODE
    import subprocess

    # get chromosome information
    cmd = "samtools idxstats ~{bamfile} | cut -f 1 | head -1"
    CHR = subprocess.check_output(cmd, shell=True)
    CHR = CHR.decode("utf-8").strip()
    print(CHR)

    coverage_file = open("~{samplename}.percent_gene_coverage.tsv", "w")
    coverage_file.write("#NOTE: THE VALUES BELOW ASSUME TBPROFILER (H37Rv) REFERENCE GENOME" + "\n")
    coverage_file.write("Gene\tPercent_Coverage" + "\n")

    with open("tbdb.bed", "r") as bedfile_fh:
      
      # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
      #  above min_depth, then wc -l counts them all
      for line in bedfile_fh:
        line = line.strip().split("\t")
        start = line[1]
        end = line[2]
        gene = line[4]

        # get depth for region in bed file 
        cmd = "samtools depth -J -r \"" + CHR + ":" + start + "-" + end + "\" ~{bamfile} | awk -F '\t' '{if (\$3 >= ~{min_depth}) print;}' | wc -l"
        print(cmd) 
        depth = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE).communicate()[0]
        print(gene + "\t" + depth)

        # get coverage for region in bed file based on depth
        # add one to gene length to compensate for subtraction
        coverage = (int(depth) / (int(end) - int(start) + 1)) * 100

        coverage_file.write(gene + "\t" + str(coverage) + "\n")
    coverage_file.close()
    CODE

    # get genome percent coverage for the entire reference genome length over min_depth
    genome=$(samtools depth -J ~{bamfile} | awk -F "\t" '{if ($3 >= ~{min_depth}) print;}' | wc -l )
    genome_pc=$(python3 -c "print ( ($genome / 4411532 ) * 100 )")
    echo "$genome_pc" > GENOME_PC
  >>>
  output {
    File tb_resistance_genes_percent_coverage = "~{samplename}.percent_gene_coverage.tsv"
    Float tb_genome_percent_coverage = read_float("GENOME_PC")
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" 
    preemptible: 0
    maxRetries: 3
  }
}