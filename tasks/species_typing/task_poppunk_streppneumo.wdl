version 1.0

task poppunk {
  meta {
    description: "Using poppunk with GPS (Global Pneumococcal Sequencing project) database for Streptococcus pneumoniae typing"
  }
  input {
    File assembly
    String samplename
    String docker = "staphb/poppunk:2.4.0"
    Int cpus = 4
    # database/reference files currently hosted on a public, requester-pays GCP bucket
    # hosting individually for speed purposes. Unzipping one big 20GB zip archive takes a long time, longer than downloading the files individually (which total 22GB uncompressed)
    # If future versions of the GPS database are released, we can update the links here or in Terra, and task should be future-proof
    File GPS_dists_npy = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.dists.npy"
    File GPS_dists_pkl = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.dists.pkl"
    File GPS_h5 = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.h5"
    File GPS_refs = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs"
    File GPS_refs_dists_npy = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.dists.npy"
    File GPS_refs_dists_pkl = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.dists.pkl"
    File GPS_refs_h5 = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6.refs.h5"
    File GPS_clusters_csv = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_clusters.csv"
    File GPS_fit_npz = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_fit.npz"
    File GPS_fit_pkl = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_fit.pkl"
    File GPS_graph_gt = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_graph.gt"
    File GPS_qcreport_txt = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_qcreport.txt"
    File GPS_unword_clusters_csv = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_unword_clusters.csv"
    File GPS_refs_graph_gt = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6refs_graph.gt"
    File GPS_external_clusters_csv = "gs://theiagen-public-files-rp/terra/theiaprok-files/GPS_v6/GPS_v6_external_clusters.csv"
  }
  command <<<
    # get version information
    poppunk --version | sed 's/poppunk //' | tee VERSION
    
    # create input TSV
    echo -e "~{samplename}\t~{assembly}" > ~{samplename}_poppunk_input.tsv
    
    # determine the database name, which is also used as a prefix for all files included in database. Also used as GPS_DB_NAME directory to put database files in
    # doing this for future proofing
    # get file name of primary h5 file, strip off suffix
    GPS_DB_NAME=$(basename ~{GPS_h5} | sed 's|.h5||')
    # sending GPS_DB_NAME into text file for logging/output purposes
    echo "${GPS_DB_NAME}" > GPS_DB_NAME

    # move all database/reference files into single directory to feed into poppunk
    mkdir -v "${GPS_DB_NAME}"
    ln -vs ~{GPS_dists_npy} ~{GPS_dists_pkl} ~{GPS_h5} ~{GPS_refs} \
      ~{GPS_refs_dists_npy} ~{GPS_refs_dists_pkl} ~{GPS_refs_h5} ~{GPS_clusters_csv} \
      ~{GPS_fit_npz} ~{GPS_fit_pkl} ~{GPS_graph_gt} ~{GPS_qcreport_txt} \
      ~{GPS_unword_clusters_csv} ~{GPS_refs_graph_gt} ~{GPS_external_clusters_csv} \
      "${GPS_DB_NAME}"/

    # to allow for compatibility with future versions of the database
    # poppunk requires this file to be explicitly passed as input
    GPS_EXTERNAL_CLUSTERS_CSV=$(ls "${GPS_DB_NAME}"/GPS_*_external_clusters.csv)

    # run poppunk
    poppunk_assign \
      --threads ~{cpus} \
      --db "${GPS_DB_NAME}" \
      --distances "${GPS_DB_NAME}/${GPS_DB_NAME}.dists" \
      --query ~{samplename}_poppunk_input.tsv \
      --output ~{samplename}_poppunk \
      --external-clustering "${GPS_EXTERNAL_CLUSTERS_CSV}"

    # parse output CSV for GPSC (Global Pneumococcal Sequence Cluster)
    if [ -f ~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv ]; then
      cut -d ',' -f 2 ~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv | tail -n 1 > GPSC.txt

      # if GPSC is "NA", overwrite with helpful message
      if [[ "$(cat GPSC.txt)" == "NA" ]]; then
        echo "Potential novel GPS Cluster identified, please email globalpneumoseq@gmail.com to have novel clusters added to the database and a GPSC cluster name assigned after you have checked for low level contamination which may contribute to biased accessory distances." >GPSC.txt
      fi
    else
      echo "poppunk failed" > GPSC.txt
    fi

  >>>
  output {
    String poppunk_gps_cluster = read_string("GPSC.txt")
    File? poppunk_gps_external_cluster_csv = "~{samplename}_poppunk/~{samplename}_poppunk_external_clusters.csv"
    String poppunk_version = read_string("VERSION")
    String poppunk_docker = docker
    String poppunk_GPS_db_version = read_string("GPS_DB_NAME")
  }
  runtime {
    docker: "~{docker}"
    # poppunk with the GPS v6 db used upwards of 12GB ram at times
    memory: "16 GB"
    cpu: cpus
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
