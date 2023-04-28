version 1.0

task snippy_variants {
  input {
    File reference_genome_file
    File read1
    File? read2
    String? query_gene
    String samplename
    String docker = "quay.io/staphb/snippy:4.6.0"
    Int cpus = 8
    Int memory = 32
    # Paramters 
    # --map_qual: Minimum read mapping quality to consider (default '60')
    # --base_quality: Minimum base quality to consider (default '13')
    # --min_coverage: Minimum site depth to for calling alleles (default '10' 
    # --min_frac: Minumum proportion for variant evidence (0=AUTO) (default '0')
    # --min_quality: Minumum QUALITY in VCF column 6 (default '100')
    # --maxsoft: Maximum soft clipping to allow (default '10')
    Int? map_qual
    Int? base_quality
    Int? min_coverage
    Float? min_frac
    Int? min_quality
    Int? maxsoft
  }
  command <<<
    snippy --version | head -1 | tee VERSION
    # set reads var
    if [ -z "~{read2}" ]; then
      reads="--se ~{read1}"
    else 
      reads="--R1 ~{read1} --R2 ~{read2}"
    fi
    # set no_hit var
    if [ -z "~{query_gene}" ]; then 
     no_hit="NA: No query gene was provided"
   else 
    no_hit="No variants identified in queried genes (~{query_gene})" 
    fi
    # call snippy
      snippy \
      --reference ~{reference_genome_file} \
      --outdir ~{samplename} \
      ${reads} \
      --cpus ~{cpus} \
      --ram ~{memory} \
      --prefix ~{samplename} \
      ~{'--mapqual ' + map_qual} \
      ~{'--basequal ' + base_quality} \
      ~{'--mincov ' + min_coverage} \
      ~{'--minfrac ' + min_frac} \
      ~{'--minqual ' + min_quality} \
      ~{'--maxsoft ' + maxsoft}
    # parse gene-specific outputs from snps.tab
    echo -e "samplename,$(head -n 1 ./~{samplename}/~{samplename}.csv)" > ./gene_query.csv
    for qgene in $(echo "~{query_gene}" | sed "s/,/ /g");
      # capture queried hits to single file 
      do 
        if grep -q  "${qgene}" ./~{samplename}/~{samplename}.csv; then 
          grep "${qgene}" ./~{samplename}/~{samplename}.csv | awk '{print "'~{samplename}'," $0}' >> ./gene_query.csv
          # curate relevant columns of quieried hits to single output
          grep "${qgene}" ./gene_query.csv | awk -F"," '{print "'${qgene}': "$15" ("$12"; "$7")"}' >> snippy_variant_hits_tmp
        fi
     done
   # convert newlines to comma

   if [ -f snippy_variant_hits_tmp ]; then
     paste -s -d, snippy_variant_hits_tmp > SNIPPY_VARIANT_HITS
   else
     echo "${no_hit}" > SNIPPY_VARIANT_HITS
   fi
    # Compress output dir
    tar -cvzf "./~{samplename}_snippy_variants_outdir.tar" "./~{samplename}"
  >>>
  output {
    String snippy_variants_version = read_string("VERSION")
    String snippy_docker = docker
    String snippy_variants_query = "~{query_gene}"
    String snippy_variants_hits = read_string("SNIPPY_VARIANT_HITS")
    File snippy_variants_outdir_tarball = "./~{samplename}_snippy_variants_outdir.tar"
    Array[File] snippy_outputs = glob("~{samplename}/~{samplename}*")
    File snippy_variants_gene_query_results = "./gene_query.csv"
    File snippy_variants_results = "~{samplename}/~{samplename}.csv"
    File snippy_variants_bam = "~{samplename}/~{samplename}.bam"
    File snippy_variants_bai ="~{samplename}/~{samplename}.bam.bai"
    File snippy_variants_summary = "~{samplename}/~{samplename}.txt"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: "~{cpus}"
      disks: "local-disk 100 SSD"
      preemptible: 0
      maxRetries: 3
  }
}