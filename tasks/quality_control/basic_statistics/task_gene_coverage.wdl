version 1.0

task gene_coverage {
  input {
    File bamfile
    File bedfile
    String samplename
    Int sc2_s_gene_start = 21563
    Int sc2_s_gene_stop = 25384
    Int min_depth = 10
    
    String organism
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/samtools:1.15"
    Int disk_size = 100
    Int memory = 8
    Int cpu = 2
  }
  command <<<
    echo "Calculating gene coverage for ~{samplename} using samtools with ~{bamfile}"

    samtools index ~{bamfile}

    # extract chromosome name from bam file
    chromosome=$(samtools idxstats ~{bamfile} | cut -f 1 | head -1)
    export chromosome

    if [ "~{organism}" == "sars-cov-2" ]; then
      echo "organism sras-"
      samtools depth -r "$chromosome:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} > ~{samplename}.s_gene.depth
      s_gene_depth=$(cut -f 7 ~{samplename}.s_gene.depth | tail -n 1)
      
      if [ -z "$s_gene_depth" ] ; then s_gene_depth="0"; fi
      echo "$s_gene_depth" | tee S_GENE_DEPTH
      
      sgene=$(samtools depth -J -r "${chromosome}:~{sc2_s_gene_start}-~{sc2_s_gene_stop}" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
      sgene_pc=$(python3 -c "print ( round( ($sgene / 3822 ) * 100, 2 ) )")
      echo "$sgene_pc" | tee S_GENE_PC
    else
      echo 0.0 | tee S_GENE_DEPTH
      echo 0.0 | tee S_GENE_PC
    fi

    # add warning to file
    echo "Caution: results may be inaccurate if your sample is not mapped to the reference genome used to generate the bed file of gene locations." > ~{samplename}.percent_gene_coverage.tsv
    # iterate through file and calculate coverage for each row in the bedfile
    while read -r line; do
      # pull out the important fields from the bedfile
      start=$(echo "$line" | cut -f 2)
      stop=$(echo "$line" | cut -f 3)
      gene=$(echo "$line" | cut -f 4)

      # calculate depth and coverage with samtools and python
      depth=$(samtools depth -r "$chromosome:$start-$stop" ~{bamfile} | awk -F "\t" '{ if ($3 > ~{min_depth}) print; }' | wc -l)
      coverage=$(python3 -c "print ( round( ($depth / ($stop - $start + 1) ) * 100, 2 ) )")

      echo -e "$gene\t$coverage" >> ~{samplename}.percent_gene_coverage.tsv
    done < ~{bedfile}

  >>>
  output {
    Float sc2_s_gene_depth = read_float("S_GENE_DEPTH")
    Float sc2_s_gene_percent_coverage = read_float("S_GENE_PC")
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