version 1.0

task usher {
  input {
    Array[File] assembly_fasta
    String tree_name
    
    # the protobuf tree used to place sequences onto -- required if not mpox or sc2
    File? mutation_annotated_tree_pb

    # the reference genome to align to -- required if not mpox or sc2
    File? reference_genome
    Int subtree_size = 20 # change value to indicate how many of the closest-related samples you want to show in a subtree

    # what organism to run usher on
    String organism # options: sars-cov-2, mpox

    # runtime
    String docker = "us-docker.pkg.dev/general-theiagen/pathogengenomics/usher:0.6.2"
    Int memory = 8
    Int cpus = 2
    Int disk_size = 100
  }
  command <<<
    if [ "~{organism}" == "mpox" ]; then
      # download latest protobuf file
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV/mpxv.latest.masked.pb.gz"
      # unzip protobuf file into renamed file
      gunzip -c mpxv.latest.masked.pb.gz > mutation.annotated.tree.pb
      
      # download versioning information of protobuf
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/UShER_hMPXV/mpxv.latest.version.txt"
      cat mpxv.latest.version.txt > PROTOBUF_VERSION

      # download mpox reference used to make usher protobuf
      wget "https://hgdownload.gi.ucsc.edu/hubs/GCF/014/621/545/GCF_014621545.1/GCF_014621545.1.fa.gz"
      # unzip reference into renamed file
      gunzip -c GCF_014621545.1.fa.gz > reference_genome.fasta
    elif [ "~{organism}" == "sars-cov-2" ]; then
      # download latest protobuf files
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
      # unzip protobuf file into renamed file
      gunzip -c public-latest.all.masked.pb.gz > mutation.annotated.tree.pb

      # download versioning information of protobuf
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt"
      cat public-latest.version.txt > PROTOBUF_VERSION

      # copy the usher SC2 reference to expected name
      cp test/NC_045512v2.fa reference_genome.fasta
    else
      # assume user provided materials
      echo "User-provided protobuf file" > PROTOBUF_VERSION
      # rename provided protobuf and reference
      mv "~{mutation_annotated_tree_pb}" mutation.annotated.tree.pb
      mv "~{reference_genome}" reference_genome.fasta
    fi

    # concatenate assembly files
    for file in ~{sep=' ' assembly_fasta}; do
      cat "${file}" >> all_sequences.fasta
    done

    # generate multiple sequence alignment
    mafft --thread 10 --auto --keeplength --addfragments all_sequences.fasta reference_genome.fasta > all_sequences.aligned.fasta

    # convert MSA to VCF
    faToVcf all_sequences.aligned.fasta all_sequences.merged.vcf

    # capture usher version
    usher --version > VERSION

    # run usher
    usher \
      --load-mutation-annotated-tree mutation.annotated.tree.pb \
      --vcf all_sequences.merged.vcf \
      --write-uncondensed-final-tree \
      --write-subtrees-size ~{subtree_size}

    mv uncondensed-final-tree.nh ~{tree_name}_uncondensed-final-tree.nwk
    mv clades.txt ~{tree_name}_clades.txt
  >>>
  output {
    File usher_uncondensed_tree = "~{tree_name}_uncondensed-final-tree.nwk"
    File usher_clades = "~{tree_name}_clades.txt"
    Array[File] usher_subtrees = glob("subtree-*.nh")
    String usher_version = read_string("VERSION")
    String usher_protobuf_version = read_string("PROTOBUF_VERSION")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu :  cpus
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x4"
  }
}