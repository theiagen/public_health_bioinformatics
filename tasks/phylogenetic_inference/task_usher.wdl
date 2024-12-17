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

    # what organism to run usher on -- the organisms below have default data provided
    String organism # options: sars-cov-2, mpox, RSV-A, RSV-B

    # runtime
    String docker = "us-docker.pkg.dev/general-theiagen/pathogengenomics/usher:0.6.2"
    Int memory = 16
    Int cpu = 4
    Int disk_size = 200
  }
  command <<<
    if [ "~{organism}" == "mpox" ]; then
      echo "DEBUG: organism is mpox, downloading latest UShER data"
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
      echo "DEBUG: finished downloading mpox data"
    elif [ "~{organism}" == "sars-cov-2" ]; then
      echo "DEBUG: organism is sars-cov-2, downloading latest UShER data"
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz"
      # unzip protobuf file into renamed file
      gunzip -c public-latest.all.masked.pb.gz > mutation.annotated.tree.pb

      # download versioning information of protobuf
      wget "https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt"
      cat public-latest.version.txt > PROTOBUF_VERSION
      
      # copy the usher SC2 reference to expected name
      cp /HOME/usher/test/NC_045512v2.fa reference_genome.fasta
      echo "DEBUG: finished downloading sars-cov-2 data"
    elif [ "~{organism}" == "RSV-A" ]; then
      echo "DEBUG: organism is RSV-A, downloading latest UShER data"
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/rsvA.latest.pb.gz"
      # unzip protobuf file into renamed file
      gunzip -c rsvA.latest.pb.gz > mutation.annotated.tree.pb

      # download versioning information of protobuf
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/UShER_RSV-A/rsvA.latest.version.txt"
      cat rsvA.latest.version.txt > PROTOBUF_VERSION

      # download RSV-A reference used to make usher protobuf
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/002/815/475/GCF_002815475.1/GCF_002815475.1.fa.gz"
      # unzip reference into renamed file
      gunzip -c GCF_002815475.1.fa.gz > reference_genome.fasta
      echo "DEBUG: finished downloading RSV-A data"
    elif [ "~{organism}" == "RSV-B" ]; then
      echo "DEBUG: organism is RSV-B, downloading latest UShER data"
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/855/545/GCF_000855545.1/UShER_RSV-B/rsvB.latest.pb.gz"
      # unzip protobuf file into renamed file
      gunzip -c rsvB.latest.pb.gz > mutation.annotated.tree.pb

      # download versioning information of protobuf
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/855/545/GCF_000855545.1/UShER_RSV-B/rsvB.latest.version.txt"
      cat rsvB.latest.version.txt > PROTOBUF_VERSION

      # download RSV-B reference used to make usher protobuf
      wget "https://hgdownload.soe.ucsc.edu/hubs/GCF/000/855/545/GCF_000855545.1/GCF_000855545.1.fa.gz"
      # unzip reference into renamed file
      gunzip -c GCF_000855545.1.fa.gz > reference_genome.fasta
      echo "DEBUG: finished downloading RSV-B data"
    else
      echo "DEBUG: organism is unknown, assuming user-provided data"
      echo "User-provided protobuf file" > PROTOBUF_VERSION
      # rename provided protobuf and reference
      mv "~{mutation_annotated_tree_pb}" mutation.annotated.tree.pb
      mv "~{reference_genome}" reference_genome.fasta
      echo "DEBUG: finished renaming user-provided data"
    fi

    echo "DEBUG: concatenating assembly files"
    for file in ~{sep=' ' assembly_fasta}; do
      cat "${file}" >> all_sequences.fasta
    done

    echo "DEBUG: generating MSA"
    mafft --thread 10 --auto --keeplength --addfragments all_sequences.fasta reference_genome.fasta > all_sequences.aligned.fasta

    echo "DEBUG: converting MSA to VCF"
    faToVcf all_sequences.aligned.fasta all_sequences.merged.vcf

    # capture usher version
    usher --version > VERSION

    echo "DEBUG: running UShER"
    usher \
      --load-mutation-annotated-tree mutation.annotated.tree.pb \
      --vcf all_sequences.merged.vcf \
      --write-uncondensed-final-tree \
      --write-subtrees-size ~{subtree_size}

    echo "DEUBG: UShER finished, renaming files"
    mv uncondensed-final-tree.nh ~{tree_name}_uncondensed-final-tree.nwk
    mv clades.txt ~{tree_name}_clades.txt

    # rename subtree files to .nwk extension
    for file in *nh; do
      mv "$file" "$(basename -- "$file" .nh).nwk"
    done

  >>>
  output {
    File usher_uncondensed_tree = "~{tree_name}_uncondensed-final-tree.nwk"
    File usher_clades = "~{tree_name}_clades.txt"
    Array[File] usher_subtrees = glob("subtree-*.nwk")
    Array[File] usher_subtree_mutations = glob("subtree-*-mutations.txt")
    String usher_version = read_string("VERSION")
    String usher_protobuf_version = read_string("PROTOBUF_VERSION")
  }
  runtime {
    docker: docker
    memory: memory + " GB"
    cpu :  cpu
    disks:  "local-disk " + disk_size + " HDD"
    disk: disk_size + " GB"
    preemptible: 0
    dx_instance_type: "mem3_ssd1_v2_x4"
  }
}
