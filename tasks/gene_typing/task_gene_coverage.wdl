version 1.0

task gene_coverage {
  input {
    File bamfile
    File baifile
    String samplename
    Int sc2_s_gene_start = 21563
    Int sc2_s_gene_stop = 25384
    Int min_depth
    
    String organism
    String docker
    Int disk_size = 100
  }
  command <<<

    if [ "~{organism}" == "sars-cov-2" ]; then
      chromosome=$(samtools idxstats ~{bamfile} | cut -f 1 | head -1)
      samtools depth -r "$chromosome:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} > ~{samplename}.s_gene.depth
      s_gene_depth=$(cut -f 7 ~{samplename}.cov.txt | tail -n 1)
      
      if [ -z "s_gene_depth" ] ; then s_gene_depth="0"; fi
      echo $s_gene_depth | tee S_GENE_DEPTH
      
      sgene=$(samtools depth -J -r "${chromosome}:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
      sgene_pc=$(python3 -c "print ( round( ($sgene / 3822 ) * 100, 2 ) )")
      echo "$sgene_pc" | tee S_GENE_PC
    fi


  >>>
  output {
    Float sc2_s_gene_depth = read_string("S_GENE_DEPTH")
    Float sc2_s_gene_percent_coverage = read_string("S_GENE_PC")
    File sc2_all_genes_percent_coverage = "~{samplename}.percent_gene_coverage.tsv"
  }
  runtime {
    docker: "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
    memory: "8 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}