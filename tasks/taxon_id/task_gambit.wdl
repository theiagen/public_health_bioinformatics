version 1.0

task gambit {
  input {
    File assembly
    String samplename
    String docker = "quay.io/staphb/gambit:0.5.0"
    File? gambit_db_genomes
    File? gambit_db_signatures
  }
  # If "File" type is used Cromwell attempts to localize it, which fails because it doesn't exist yet.
  String report_path = "~{samplename}_gambit.json"
  String closest_genomes_path = "~{samplename}_gambit_closest.csv"
  command <<<
    # capture date and version
    date | tee DATE
    gambit --version | tee GAMBIT_VERSION
    
    # set gambit reference dir; will assume that gambit genomes and signatures will be provided by user in tandem or not at all
    if [[ ! -z "~{gambit_db_genomes}" ]]; then 
      echo "User gabmit db identified; ~{gambit_db_genomes} will be utilized for alignment"
      gambit_db_version="$(basename -- '~{gambit_db_genomes}'); $(basename -- '~{gambit_db_signatures}')"
      gambit_db_dir="${PWD}/gambit_database"
      mkdir ${gambit_db_dir}
      cp ~{gambit_db_genomes} ${gambit_db_dir}
      cp ~{gambit_db_signatures} ${gambit_db_dir}
    else
     gambit_db_dir="/gambit-db" 
     gambit_db_version="unmodified from gambit container: ~{docker}"
    fi
    
    echo ${gambit_db_version} | tee GAMBIT_DB_VERSION
    
    gambit -d ${gambit_db_dir} query -f json -o ~{report_path} ~{assembly} 
    
    python3 <<EOF
    import json
    import csv

    def fmt_dist(d): return format(d, '.4f')

    with open("~{report_path}") as f:
      data = json.load(f)

    (item,) = data['items']
    predicted = item['predicted_taxon']
    next_taxon = item['next_taxon']
    closest = item['closest_genomes'][0]

    with open('CLOSEST_DISTANCE', 'w') as f:
      f.write(fmt_dist(closest['distance']))

    # Predicted taxon
    with open('PREDICTED_TAXON', 'w') as f:
      if predicted is None:
        f.write('NA')
      elif predicted['name'] is None:
        f.write('NA')
      else:
        f.write(predicted['name'])
    with open('PREDICTED_TAXON_RANK', 'w') as f:
      if predicted is None:
        f.write('NA')
      elif predicted['rank'] is None:
        f.write('NA')
      else:
        f.write(predicted['rank'])
    with open('PREDICTED_TAXON_THRESHOLD', 'w') as f:
      if predicted is None:
        f.write(fmt_dist(0))
      elif predicted['distance_threshold'] is None:
        f.write(fmt_dist(0))
      else:
        f.write(fmt_dist(predicted['distance_threshold']))

    # Next taxon
    with open('NEXT_TAXON', 'w') as f:
      if next_taxon is None:
        f.write('NA')
      elif next_taxon['name'] is None:
        f.write('NA')
      else:
        f.write(next_taxon['name'])
    with open('NEXT_TAXON_RANK', 'w') as f:
      if next_taxon is None:
        f.write('NA')
      elif next_taxon['rank'] is None:
        f.write('NA')
      else:
        f.write(next_taxon['rank'])
    with open('NEXT_TAXON_THRESHOLD', 'w') as f:
      if next_taxon is None:
        f.write(fmt_dist(0))
      elif next_taxon['distance_threshold'] is None:
        f.write(fmt_dist(0))
      else:
        f.write(fmt_dist(next_taxon['distance_threshold']))
      
    # Table of closest genomes
    with open('~{closest_genomes_path}', 'w', newline='') as f:
      writer = csv.writer(f)

      # Header
      writer.writerow([
        'distance',
        'genome.description',
        'genome.taxon.name',
        'genome.taxon.rank',
        'genome.taxon.threshold',
        'matched.name',
        'matched.rank',
        'matched.distance_threshold',
      ])

      for match in item['closest_genomes']:
        genome = match['genome']
        genome_taxon = genome['taxonomy'][0]
        match_taxon = match['matched_taxon']

        writer.writerow([
          fmt_dist(match['distance']),
          genome['description'],
          genome_taxon['name'],
          genome_taxon['rank'],
          fmt_dist(genome_taxon['distance_threshold']),
          '' if match_taxon is None else match_taxon['name'],
          '' if match_taxon is None else match_taxon['rank'],
          fmt_dist(0 if match_taxon is None else match_taxon['distance_threshold']),
        ])
    EOF
    # set merlin tags
    predicted_taxon=$(cat PREDICTED_TAXON)
    if [[ ${predicted_taxon} == *"Escherichia"* ]] || [[ ${predicted_taxon} == *"Shigella"* ]] ; then 
      merlin_tag="Escherichia"
    elif [[ ${predicted_taxon} == *"Haemophilus"* ]]; then 
      merlin_tag="Haemophilus"
    elif [[ ${predicted_taxon} == *"Klebsiella"* ]]; then 
      merlin_tag="Klebsiella"
    elif [[ ${predicted_taxon} == *"Acinetobacter baumannii"* ]]; then 
      merlin_tag="Acinetobacter baumannii"
    elif [[ ${predicted_taxon} == *"Legionella pneumophila"* ]]; then 
      merlin_tag="Legionella pneumophila"
    elif [[ ${predicted_taxon} == *"Listeria"* ]]; then 
      merlin_tag="Listeria"
    elif [[ ${predicted_taxon} == *"Mycobacterium tuberculosis"* ]]; then 
      merlin_tag="Mycobacterium tuberculosis"
    elif [[ ${predicted_taxon} == *"Neisseria"* ]]; then 
      merlin_tag="Neisseria"
    elif [[ ${predicted_taxon} == *"Salmonella"* ]]; then 
      merlin_tag="Salmonella"
    elif [[ ${predicted_taxon} == *"Staphylococcus"* ]]; then 
      merlin_tag="Staphylococcus"
    elif [[ ${predicted_taxon} == *"Streptococcus"* ]]; then 
      merlin_tag="Streptococcus"
      # set to pneumoniae if gambit calls the species
      if [[ ${predicted_taxon} == *"Streptococcus pneumoniae"* ]]; then 
        merlin_tag="Streptococcus pneumoniae"
      fi
    else 
      merlin_tag="None"
    fi
    echo ${merlin_tag} | tee MERLIN_TAG
  >>>
  output {
    File gambit_report_file = report_path
    File gambit_closest_genomes_file = closest_genomes_path
    String gambit_predicted_taxon = read_string("PREDICTED_TAXON")
    String gambit_predicted_taxon_rank = read_string("PREDICTED_TAXON_RANK") 
    String gambit_next_taxon = read_string("NEXT_TAXON")
    String gambit_next_taxon_rank = read_string("NEXT_TAXON_RANK")
    String gambit_version = read_string("GAMBIT_VERSION")
    String gambit_db_version = read_string("GAMBIT_DB_VERSION")
    String merlin_tag = read_string("MERLIN_TAG")
    String gambit_docker = docker
  }
  runtime {
    docker:  "~{docker}"
    memory:  "16 GB"
    cpu:   8
    disks: "local-disk 100 SSD"
    preemptible:  0
  }
}
