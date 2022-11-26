version 1.0

task shovill_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "quay.io/staphb/shovill:1.1.0"

    ## SHOVILL optional parameters
    ##  --depth [INT]           Sub-sample --R1/--R2 to this depth. Disable with --depth 0 (default: 150)
    ##  --gsize [STRING]        Estimated genome size eg. 3.2M <blank=AUTODETECT> (default: '')
    ##  --minlen [INT]          Minimum contig length <0=AUTO> (default: 0)
    ##  --mincov [FLOAT]        Minimum contig coverage <0=AUTO> (default: 2)
    ##  --assembler [STRING]    Assembler: skesa velvet megahit spades (default: 'spades')
    ##  --opts [STRING]         Extra assembler options in quotes eg. spades: "--untrusted-contigs locus.fna" ... (default: '')
    ##  --kmers [STRING]        K-mers to use <blank=AUTO> (default: '')
    ##  --trim [BOOLEAN]        Enable adaptor trimming (default: OFF)
    ##  --noreadcorr [BOOLEAN]  Disable read error correction (default: OFF)
    ##  --nostitch [BOOLEAN]    Disable read stitching (default: OFF)
    ##  --nocorr [BOOLEAN]      Disable post-assembly correction (default: OFF)

    
    Int? depth
    String? genome_size
    Int min_contig_length = 200
    Float? min_coverage
    String assembler = "skesa"
    String? assembler_options
    String? kmers
    Boolean trim = false
    Boolean noreadcorr = false
    Boolean nostitch = false
    Boolean nocorr = false
  }
  command <<<
    shovill --version | head -1 | tee VERSION
    shovill \
    --outdir out \
    --R1 ~{read1_cleaned} \
    --R2 ~{read2_cleaned} \
    --minlen ~{min_contig_length} \
    ~{'--depth ' + depth} \
    ~{'--gsize ' + genome_size} \
    ~{'--mincov ' + min_coverage} \
    ~{'--assembler ' + assembler} \
    ~{'--opts ' + assembler_options} \
    ~{'--kmers ' + kmers} \
    ~{true='--trim' false='' trim} \
    ~{true='--noreadcorr' false='' noreadcorr} \
    ~{true='--nostitch' false='' nostitch} \
    ~{true='--nocorr' false='' nocorr}

    mv out/contigs.fa out/~{samplename}_contigs.fasta

    if [ "~{assembler}" == "spades" ] ; then
      mv out/contigs.gfa out/~{samplename}_contigs.gfa
    elif [ "~{assembler}" == "megahit" ] ; then
      mv out/contigs.fastg out/~{samplename}_contigs.fastg
    elif [ "~{assembler}" == "velvet" ] ; then
      mv out/contigs.LastGraph out/~{samplename}_contigs.LastGraph
    fi
    
  >>>
  output {
    File assembly_fasta = "out/~{samplename}_contigs.fasta"
    File? contigs_gfa = "out/~{samplename}_contigs.gfa"
    File? contigs_fastg = "out/~{samplename}_contigs.fastg"
    File? contigs_lastgraph = "out/~{samplename}_contigs.LastGraph"
    String shovill_version = read_string("VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "16 GB"
      cpu: 4
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}

task shovill_se {
  input {
    File read1_cleaned
    String samplename
    String docker = "quay.io/staphb/shovill-se:1.1.0"

    ## SHOVILL optional parameters
    ##  --depth [INT]           Sub-sample --R1/--R2 to this depth. Disable with --depth 0 (default: 150)
    ##  --gsize [STRING]        Estimated genome size eg. 3.2M <blank=AUTODETECT> (default: '')
    ##  --minlen [INT]          Minimum contig length <0=AUTO> (default: 0)
    ##  --mincov [FLOAT]        Minimum contig coverage <0=AUTO> (default: 2)
    ##  --assembler [STRING]    Assembler: skesa velvet megahit spades (default: 'spades')
    ##  --opts [STRING]         Extra assembler options in quotes eg. spades: "--untrusted-contigs locus.fna" ... (default: '')
    ##  --kmers [STRING]        K-mers to use <blank=AUTO> (default: '')
    ##  --trim [BOOLEAN]        Enable adaptor trimming (default: OFF)
    ##  --noreadcorr [BOOLEAN]  Disable read error correction (default: OFF)
    ##  --nocorr [BOOLEAN]      Disable post-assembly correction (default: OFF)

    Int? depth
    String? genome_size
    Int min_contig_length = 200
    Float? min_coverage
    String assembler = "spades"
    String? assembler_options
    String? kmers
    Boolean trim = false
    Boolean noreadcorr = false
    Boolean nocorr = false
  }
  command <<<
    shovill-se --version | head -1 | tee VERSION
    shovill-se \
    --outdir out \
    --se ~{read1_cleaned} 
    --minlen ~{min_contig_length} \
    ~{'--depth ' + depth} \
    ~{'--gsize ' + genome_size} \
    ~{'--mincov ' + min_coverage} \
    ~{'--assembler ' + assembler} \
    ~{'--opts ' + assembler_options} \
    ~{'--kmers ' + kmers} \
    ~{true='--trim' false='' trim} \
    ~{true='--noreadcorr' false='' noreadcorr} \
    ~{true='--nocorr' false='' nocorr}

    mv out/contigs.fa out/~{samplename}_contigs.fasta

    if [ "~{assembler}" == "spades" ] ; then
      mv out/contigs.gfa out/~{samplename}_contigs.gfa
    elif [ "~{assembler}" == "megahit" ] ; then
      mv out/contigs.fastg out/~{samplename}_contigs.fastg
    elif [ "~{assembler}" == "velvet" ] ; then
      mv out/contigs.LastGraph out/~{samplename}_contigs.LastGraph
    fi
  >>>
  output {
    File assembly_fasta = "out/~{samplename}_contigs.fasta"
    File? contigs_gfa = "out/~{samplename}_contigs.gfa"
    File? contigs_fastg = "out/~{samplename}_contigs.fastg"
    File? contigs_lastgraph = "out/~{samplename}_contigs.LastGraph"
    String shovill_version = read_string("VERSION")
  }
  runtime {
      docker: "~{docker}"
      memory: "16 GB"
      cpu: 4
      disks: "local-disk 100 SSD"
      preemptible: 0
  }
}
