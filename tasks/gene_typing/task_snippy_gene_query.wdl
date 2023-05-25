version 1.0

task snippy_gene_query {
  input {
    String samplename
    File snippy_variants_results
    String? query_gene
    File? reference
    String docker = "quay.io/staphb/snippy:4.6.0"
    Int cpus = 8
    Int memory = 32
  }
  command <<<
    # set variable for if no hits are detected for query genes
    if [ -z "~{query_gene}" ]; then 
        no_hit="NA: No query gene was provided"
    else 
        no_hit="No variants identified in queried genes (~{query_gene})" 
    fi

    # if provided, check that query gene strings are present in reference genome
    if [ -z "~{reference}" ]; then 
      echo "No reference genome was provided, query genes were not verified to be in reference genome" >> QUERY_CHECK
    else 
      echo "DEBUG: Reference genome was provided, checking that gene queries are in reference genome" 
      for qgene in $(echo "~{query_gene}" | sed "s/,/ /g"); do 
          echo "DEBUG: checking reference genome for ${qgene}"
          num_lines=$(grep "${qgene}" ~{reference} | wc -l)
          if [ $num_lines -gt 0 ]; then               
            echo "${QUERY_CHECK}${qgene} was found in reference genome" >> QUERY_CHECK
          else
            echo "${QUERY_CHECK}${qgene} was NOT found in reference genome" >> QUERY_CHECK
          fi
      done
    fi

    paste -s -d, QUERY_CHECK > QUERY_CHECK_RESULTS

    # parse gene-specific outputs from snps.tab
    echo -e "samplename,$(head -n 1 ~{snippy_variants_results})" > ./gene_query.csv
    for qgene in $(echo "~{query_gene}" | sed "s/,/ /g");
      # capture queried hits to single file 
      do 
        if grep -q  "${qgene}" ~{snippy_variants_results}; then 
          grep "${qgene}" ~{snippy_variants_results} | awk '{print "'~{samplename}'," $0}' >> ./gene_query.csv
          # curate relevant columns of queried hits to single output
          grep "${qgene}" ./gene_query.csv | awk -F"," '{print "'${qgene}': "$15" ("$12"; "$7")"}' >> snippy_variant_hits_tmp
        fi
     done

    # convert newlines to comma
    if [ -f snippy_variant_hits_tmp ]; then
      paste -s -d, snippy_variant_hits_tmp > SNIPPY_VARIANT_HITS
    else
      echo "${no_hit}" > SNIPPY_VARIANT_HITS
    fi

  >>>
  output {
    String snippy_variants_query = "~{query_gene}"
    String snippy_variants_hits = read_string("SNIPPY_VARIANT_HITS")
    String snippy_variants_query_check = read_string("QUERY_CHECK_RESULTS")
    File snippy_variants_gene_query_results = "./gene_query.csv"
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