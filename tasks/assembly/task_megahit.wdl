version 1.0

task megahit_pe {
  input {
    File read1_cleaned
    File read2_cleaned
    String samplename
    String docker = "quay.io/biocontainers/megahit:1.2.9--h43eeafb_4"
    Int disk_size = 100
    Int cpu = 4
    Int memory = 16

    ## MEGAHIT optional parameters
    ##  --min-count [INT]          minimum multiplicity for filtering (k_min+1)-mers [2]
    ##  --k-list [INT,INT,..]      comma-separated list of kmer size
    ##                             all must be odd, in the range 15-255, increment <= 28)
    ##                             [21,29,39,59,79,99,119,141]
    ##  --no-mercy []              do not add mercy kmers
    ##  --bubble-level [INT]       intensity of bubble merging (0-2), 0 to disable [2]
    ##  --merge-level [L,S]        merge complex bubbles of length <= l*kmer_size and similarity >= s [20,0.95]
    ##  --prune-level [INT]        strength of low depth pruning (0-3) [2]
    ##  --prune-depth [INT]        remove unitigs with avg kmer depth less than this value [2]
    ##  --disconnect-ratio [FLOAT] disconnect unitigs if its depth is less than this ratio times
    ##                             the total depth of itself and its siblings [0.1]
    ##  --low-local-ratio [FLOAT]  remove unitigs if its depth is less than this ratio times
    ##                             the average depth of the neighborhoods [0.2]
    ##  --max-tip-len [INT]        remove tips less than this value [2*k]
    ##  --cleaning-rounds [INT]    number of rounds for graph cleanning [5]
    ##  --no-local []              disable local assembly
    ##  --kmin-1pass []            use 1pass mode to build SdBG of k_min
    ##  --presets [STRING]         override a group of parameters; possible values:
    ##                             meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'
    ##                             meta-large: '--k-min 27 --k-max 127 --k-step 10'
    ##                             (large & complex metagenomes, like soil)
    ##  -m/--memory [FLOAT]        max memory in byte to be used in SdBG construction
    ##                             (if set between 0-1, fraction of the machine's total memory) [0.9]
    ##  --mem-flag [INT]           SdBG builder memory mode. 0: minimum; 1: moderate;
    ##                             others: use all memory specified by '-m/--memory' [1]
    ##  -t/--num-cpu-threads [INT] number of CPU threads [# of logical processors]
    ##  --no-hw-accel []           run MEGAHIT without BMI2 and POPCNT hardware instructions
    ##  -o/--out-dir [STRING]      output directory [./megahit_out]  
    ##  --out-prefix [STRING]      output prefix (the contig file will be OUT_DIR/OUT_PREFIX.contigs.fa)
    ##  --min-contig-len [INT]     minimum length of contigs to output [200]
    ##  --keep-tmp-files []        keep all temporary files
    ##  --tmp-dir [STRING]         set temp directory
    
    Int min_contig_length = 200
    String? kmers
    Boolean keep_temp_files = false
  }
  command <<<
    megahit --version | head -1 | cut -d ' ' -f 2 | tee VERSION
    megahit \
      -1 ~{read1_cleaned} \
      -2 ~{read2_cleaned} \
      --min-contig-len ~{min_contig_length} \
      ~{'--k-list ' + kmers} \
      ~{true='--keep-tmp-files' false='' keep_temp_files} \
      -m ~{memory}
      -t ~{cpu}

    mv megahit_out/final.contigs.fa ~{samplename}_contigs.fasta

  >>>
  output {
    File assembly_fasta = "~{samplename}_contigs.fasta"
    String shovill_version = read_string("VERSION")
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 0
  }
}

