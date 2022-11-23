version 1.0

task snippy_variants {
  input {
    File reference
    File read1
    File? read2
    String? query_gene
    String samplename
    String docker = "staphb/snippy:4.6.0"
    Int cpus = 16
    Int memory = 64
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
    # call snippy
      snippy \
      --reference ~{reference} \
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
    echo -e "sample\t$(head -n 1 ./~{samplename}/~{samplename}.tab)" > ./gene_query.tsv
    for qgene in $(echo "~{query_gene}" | sed "s/,/ /g");
      # capture queried hits to single file 
      do 
        if grep -q  "${qgene}" ./~{samplename}/~{samplename}.tab; then 
          echo -e "samplename\t$(grep ${qgene} ./~{samplename}/~{samplename}.tab)" >> ./gene_query.tsv
          # curate relevant columns of quieried hits to single output
          grep "${qgene}" ./gene_query.tsv | awk -F"\t" '{print "'${qgene}': "$15" ("$12"; "$7")"}' >> snippy_variant_hits_tmp
        fi
     done
   # convert newlines to comma
   if [ -f snippy_variant_hits_tmp ]; then
     paste -s -d, snippy_variant_hits_tmp > SNIPPY_VARIANT_HITS
   else
     echo "No variants identified in quieried genes (~{query_gene})" > SNIPPY_VARIANT_HITS
   fi
  >>>
  output {
    String snippy_variants_version = read_string("VERSION")
    String snippy_variants_query = "~{query_gene}"
    String snippy_variants_hits = read_string("SNIPPY_VARIANT_HITS")
    File snippy_variants_gene_query_results = "./gene_query.tsv"
    Array[File] snippy_outputs =  glob("~{samplename}/~{samplename}*")
    File snippy_variants_results = "~{samplename}/~{samplename}.csv"
  }
  runtime {
      docker: "~{docker}"
      memory: "~{memory} GB"
      cpu: "~{cpus}"
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}