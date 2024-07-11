version 1.0

task genoflu {
  input {
    File assembly_fasta
    String samplename

    # excel file to cross-reference BLAST findings; probably useful if novel
    #  genotypes are not in the default file used by genoflu.py
    File? cross_reference

    Int cpu = 1
    Int disk_size = 25
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/genoflu:1.03"
    Int memory = 2
  }
  command <<<
   
    cp ~{assembly_fasta} .

    genoflu.py -v | sed -e 's/genoflu.py:\ version\ //' > VERSION

    genoflu.py \
      --fasta ~{assembly_fasta} \
      --sample_name ~{samplename} \
      ~{"--cross_reference" + cross_reference} > genoflu.output.txt

    grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2 | cut -d " " -f2 | cut -d ":" -f1 > GENOTYPE
    grep "~{samplename} Genotype" genoflu.output.txt | cut -d ">" -f2 | cut -d " " -f3- > ALL_SEGMENTS
    
    mv ~{samplename}_*_stats.tsv ~{samplename}_stats.tsv
  >>>
  output {
    String genoflu_version = read_string("VERSION")
    String genoflu_genotype = read_string("GENOTYPE")
    String genoflu_all_segments = read_string("ALL_SEGMENTS")
    File genoflu_output_tsv = "~{samplename}_stats.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: memory + " GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}