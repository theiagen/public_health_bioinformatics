version 1.0

task snippy_variants {
  input {
    File reference_genome_file
    File read1
    File? read2
    String samplename
    String docker = "quay.io/staphb/snippy:4.6.0"
    Int cpus = 8
    Int memory = 32
    # Paramters 
    # --map_qual: Minimum read mapping quality to consider (default '60')
    # --base_quality: Minimum base quality to consider (default '13')
    # --min_coverage: Minimum site depth to for calling alleles (default '10') 
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

    # Compress output dir
    tar -cvzf "./~{samplename}_snippy_variants_outdir.tar" "./~{samplename}"
  >>>
  output {
    String snippy_variants_version = read_string("VERSION")
    String snippy_variants_docker = docker
    File snippy_variants_outdir_tarball = "./~{samplename}_snippy_variants_outdir.tar"
    Array[File] snippy_variants_outputs = glob("~{samplename}/~{samplename}*")
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