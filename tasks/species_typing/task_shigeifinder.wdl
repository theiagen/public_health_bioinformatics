version 1.0

task shigeifinder {
  input {
    File assembly
    String samplename
    String docker = "staphb/shigeifinder:1.3.3"
    Int disk_size = 100
    Int cpu = 2
  }
  command <<<
    # capture date
    date | tee DATE
    # shigeifinder does not have a --version flag, relying upon the docker image tag for the version - which StaPH-B/Curtis maintains
    echo "~{docker}" | sed 's|staphb/shigeifinder:||' | tee VERSION.txt

    # ShigEiFinder checks that all dependencies are installed before running
    echo "checking for shigeifinder dependencies and running ShigEiFinder..."
    # run shigeifinder on assembly; default is 4cpus, so turning down to 2 since it's already very fast
    shigeifinder -i ~{assembly} \
        -t ~{cpu} \
        --hits \
        --output ~{samplename}_shigeifinder.tsv

    # parse output TSV
    echo "Parsing ShigEiFinder output TSV..."
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 2 >shigeifinder_ipaH_presence_absence.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 3 >shigeifinder_num_virulence_plasmid_genes.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 4 >shigeifinder_cluster.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 5 >shigeifinder_serotype.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 6 >shigeifinder_O_antigen.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 7 >shigeifinder_H_antigen.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 8 >shigeifinder_notes.txt

    # set helpful output strings if field in TSV is blank by overwriting output TXT files
    if [ "$(cat shigeifinder_ipaH_presence_absence.txt)" == "" ]; then
       echo "ShigEiFinder ipaH field was empty" > shigeifinder_ipaH_presence_absence.txt
    fi
    if [ "$(cat shigeifinder_num_virulence_plasmid_genes.txt)" == "" ]; then
       echo "ShigEiFinder number of virulence plasmid genes field was empty" > shigeifinder_num_virulence_plasmid_genes.txt
    fi
    if [ "$(cat shigeifinder_cluster.txt)" == "" ]; then
       echo "ShigEiFinder cluster field was empty" > shigeifinder_cluster.txt
    fi
    if [ "$(cat shigeifinder_serotype.txt)" == "" ]; then
       echo "ShigEiFinder serotype field was empty" > shigeifinder_serotype.txt
    fi
    if [ "$(cat shigeifinder_O_antigen.txt)" == "" ]; then
       echo "ShigEiFinder O antigen field was empty" > shigeifinder_O_antigen.txt
    fi
    if [ "$(cat shigeifinder_H_antigen.txt)" == "" ]; then
       echo "ShigEiFinder H antigen field was empty" > shigeifinder_H_antigen.txt
    fi
    if [ "$(cat shigeifinder_notes.txt)" == "" ]; then
       echo "ShigEiFinder notes field was empty" > shigeifinder_notes.txt
    fi

  >>>
  output {
    File shigeifinder_report = "~{samplename}_shigeifinder.tsv"
    String shigeifinder_docker = docker
    String shigeifinder_version = read_string("VERSION.txt")
    String shigeifinder_ipaH_presence_absence = read_string("shigeifinder_ipaH_presence_absence.txt")
    String shigeifinder_num_virulence_plasmid_genes = read_string("shigeifinder_num_virulence_plasmid_genes.txt")
    String shigeifinder_cluster = read_string("shigeifinder_cluster.txt")
    String shigeifinder_serotype = read_string("shigeifinder_serotype.txt")
    String shigeifinder_O_antigen = read_string("shigeifinder_O_antigen.txt")
    String shigeifinder_H_antigen = read_string("shigeifinder_H_antigen.txt")
    String shigeifinder_notes = read_string("shigeifinder_notes.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}
task shigeifinder_reads {
  input {
    File read1
    File? read2
    String samplename
    String docker = "staphb/shigeifinder:1.3.3"
    Int disk_size = 100
    Int cpu = 4
    Boolean paired_end = true
  }
  command <<<
    # capture date
    date | tee DATE
    # shigeifinder does not have a --version flag, relying upon the docker image tag for the version - which StaPH-B/Curtis maintains
    echo "~{docker}" | sed 's|staphb/shigeifinder:||' | tee VERSION.txt

    # ShigEiFinder checks that all dependencies are installed before running
    echo "checking for shigeifinder dependencies and running ShigEiFinder..."
    # run shigeifinder on reads; default is 4cpus, so keeping at 4 since it's doing alignment
    shigeifinder -r -i ~{read1} ~{read2} \
        ~{true='' false='--single_end' paired_end} \
        -t ~{cpu} \
        --hits \
        --output ~{samplename}_shigeifinder.tsv

    # parse output TSV
    echo "Parsing ShigEiFinder output TSV..."
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 2 >shigeifinder_ipaH_presence_absence.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 3 >shigeifinder_num_virulence_plasmid_genes.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 4 >shigeifinder_cluster.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 5 >shigeifinder_serotype.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 6 >shigeifinder_O_antigen.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 7 >shigeifinder_H_antigen.txt
    head -n 2 ~{samplename}_shigeifinder.tsv | tail -n 1 | cut -f 8 >shigeifinder_notes.txt

    # set helpful output strings if field in TSV is blank by overwriting output TXT files
    if [ "$(cat shigeifinder_ipaH_presence_absence.txt)" == "" ]; then
       echo "ShigEiFinder ipaH field was empty" > shigeifinder_ipaH_presence_absence.txt
    fi
    if [ "$(cat shigeifinder_num_virulence_plasmid_genes.txt)" == "" ]; then
       echo "ShigEiFinder number of virulence plasmid genes field was empty" > shigeifinder_num_virulence_plasmid_genes.txt
    fi
    if [ "$(cat shigeifinder_cluster.txt)" == "" ]; then
       echo "ShigEiFinder cluster field was empty" > shigeifinder_cluster.txt
    fi
    if [ "$(cat shigeifinder_serotype.txt)" == "" ]; then
       echo "ShigEiFinder serotype field was empty" > shigeifinder_serotype.txt
    fi
    if [ "$(cat shigeifinder_O_antigen.txt)" == "" ]; then
       echo "ShigEiFinder O antigen field was empty" > shigeifinder_O_antigen.txt
    fi
    if [ "$(cat shigeifinder_H_antigen.txt)" == "" ]; then
       echo "ShigEiFinder H antigen field was empty" > shigeifinder_H_antigen.txt
    fi
    if [ "$(cat shigeifinder_notes.txt)" == "" ]; then
       echo "ShigEiFinder notes field was empty" > shigeifinder_notes.txt
    fi

  >>>
  output {
    File shigeifinder_report = "~{samplename}_shigeifinder.tsv"
    String shigeifinder_docker = docker
    String shigeifinder_version = read_string("VERSION.txt")
    String shigeifinder_ipaH_presence_absence = read_string("shigeifinder_ipaH_presence_absence.txt")
    String shigeifinder_num_virulence_plasmid_genes = read_string("shigeifinder_num_virulence_plasmid_genes.txt")
    String shigeifinder_cluster = read_string("shigeifinder_cluster.txt")
    String shigeifinder_serotype = read_string("shigeifinder_serotype.txt")
    String shigeifinder_O_antigen = read_string("shigeifinder_O_antigen.txt")
    String shigeifinder_H_antigen = read_string("shigeifinder_H_antigen.txt")
    String shigeifinder_notes = read_string("shigeifinder_notes.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
    maxRetries: 3
  }
}