version 1.0

task gene_coverage {
  input {
    File bamfile
    File baifile
    File bedfile
    String samplename
    Int sc2_s_gene_start = 21563
    Int sc2_s_gene_stop = 25384
    Int min_depth
    
    String organism
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    echo "Calculating gene coverage for ~{samplename} using samtools with ~{bamfile} and ~{baifile}"

    if [ "~{organism}" == "sars-cov-2" ]; then
      chromosome=$(samtools idxstats ~{bamfile} | cut -f 1 | head -1)
      samtools depth -r "$chromosome:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} > ~{samplename}.s_gene.depth
      s_gene_depth=$(cut -f 7 ~{samplename}.cov.txt | tail -n 1)
      
      if [ -z "s_gene_depth" ] ; then s_gene_depth="0"; fi
      echo $s_gene_depth | tee S_GENE_DEPTH
      
      sgene=$(samtools depth -J -r "${chromosome}:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
      sgene_pc=$(python3 -c "print ( round( ($sgene / 3822 ) * 100, 2 ) )")
      echo "$sgene_pc" | tee S_GENE_PC
    else
      echo "0" | tee S_GENE_DEPTH
      echo "0" | tee S_GENE_PC
    fi

    python3 <<CODE
      import subprocess

      with open(~{bedfile}, "r") as bedfile, open("~{samplename}.percent_gene_coverage.tsv", "w") as outfile:
        outfile.write("Gene\Percent_Coverage\n")
        for line in bedfile:
          if line.startswith("#"): continue
          line = line.strip().split("\t")
          chromosome = line[0]
          start = line[1]
          stop = line[2]
          gene = line[3]

          command = "samtools depth -r \"" + chromosome + ":" + start + "-" + stop + "\" " + {bamfile} + " | wc -l "
          depth = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=True).communicate()[0]
        
          # get coverage
          coverage = (int(depth) / (int(stop) - int(start) + 1)) * 100
          outfile.write("{}\t{}\n".format(gene, coverage))

    CODE

  >>>
  output {
    Float sc2_s_gene_depth = read_string("S_GENE_DEPTH")
    Float sc2_s_gene_percent_coverage = read_string("S_GENE_PC")
    File est_percent_gene_coverage_tsv = "~{samplename}.percent_gene_coverage.tsv"
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}