version 1.0

task usher {
  input {
    Array[File] assembly_fasta
    String tree_name
    
    # the protobuf tree used to place sequences onto
    File? mutation_annotated_tree_pb

    # the reference genome to align to
    File? reference_genome

    # what organism to run usher on
    String organism = "mpox" # options: sars-cov-2, mpox, ??

    # runtime
    String docker = "us-docker.pkg.dev/general-theiagen/pathogengenomics/usher:0.6.2"
    Int memory = 8
    Int cpus = 2
    Int disk_size = 100
  }
  command <<<
    if [ "~{organism}" == "mpox" ]; then
      # download latest protobuf files
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV/mpxv.latest.masked.pb.gz"
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV/mpxv.latest.version.txt"

      # download mpox reference used to make usher protobuf
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/GCF_014621545.1.fa.gz"
      # unzip reference
      gunzip GCF_014621545.1.fa.gz 
      # rename reference file
      mv GCF_014621545.1.fa reference_genome.fasta

      # capture version of protobuf
      cat mpxv.latest.version.txt > PROTOBUF_VERSION

      # unzip protobuf file
      gunzip mpxv.latest.masked.pb.gz
      # rename protobuf 
      mv mpxv.latest.masked.pb mutation.annotated.tree.pb
    elif [ "~{organism}" == "sars-cov-2" ]; then
      # download latest protobuf files
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt"

      # copy the usher SC2 reference to expected name
      cp test/NC_045512v2.fa reference_genome.fasta

      # capture version of protobuf
      cat public-latest.version.txt > PROTOBUF_VERSION

      # unzip protobuf file
      gunzip public-latest.all.masked.pb.gz
      # rename protobuf
      mv public-latest.all.masked.pb mutation.annotated.tree.pb
    else
      # assume user provided materials
      echo "User-provided protobuf file" > PROTOBUF_VERSION
      # rename provided protobuf and reference
      mv "~{mutation_annotated_tree_pb}" mutation.annotated.tree.pb
      mv "~{reference_genome}" reference_genome.fasta
    fi

    # concatenate assembly files
    for file in assembly_fasta; do
      cat $file >> all_sequences.fasta
    done

    # generate multiple sequence alignment
    mafft --thread 10 --auto --keeplength --addfragments all_sequences.fasta reference_genome.fasta > all_sequences.aligned.fasta

    # convert MSA to VCF
    faToVcf <(cat reference_genome.fasta all_sequences.aligned.fasta) all_sequences.merged.vcf

    # run usher
    usher -i mutation.annotated.tree.pb -v all_sequences.merged.vcf -u 

    mv uncondensed-final-tree.nh ~{tree_name}_final-tree.nwk
  >>>
  output {
    File usher_tree = "~{tree_name}_final-tree.nwk"
    String usher_version = read_string("VERSION")
    String mpox_protobuf_version = read_string("PROTOBUF_VERSION")

  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu :  cpus
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x4"
    maxRetries: 3
  }
}