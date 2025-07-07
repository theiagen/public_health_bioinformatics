version 1.0

task gambit {
  input {
    File assembly
    String samplename
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0"
    File gambit_db_genomes = "gs://gambit-databases-rp/2.0.0/gambit-metadata-2.0.1-20250505.gdb"
    File gambit_db_signatures = "gs://gambit-databases-rp/2.0.0/gambit-signatures-2.0.1-20250505.gs"
    Int disk_size = 20
    Int memory = 2
    Int cpu = 1
  }
  # If "File" type is used Cromwell attempts to localize it, which fails because it doesn't exist yet.
  String report_path = "~{samplename}_gambit.json"
  String closest_genomes_path = "~{samplename}_gambit_closest.csv"
  command <<<
    # capture date and version
    date | tee DATE
    gambit --version | tee GAMBIT_VERSION
    
    # set gambit reference dir; will assume that gambit genomes and signatures will be provided by user in tandem or not at all
    # -s evaluates to TRUE if the file exists and has a size greater than zero
    if [[ -s "~{gambit_db_genomes}" ]]; then 
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
    
    gambit -d ${gambit_db_dir} query -f json -o ~{report_path} ~{assembly} -c ~{cpu}
    
    python3 <<EOF
    import json
    import csv
    import re

    def fmt_dist(d): return format(d, '.4f')

    with open("~{report_path}") as f:
      data = json.load(f)

    (item,) = data['items']
    predicted = item['predicted_taxon']
    next_taxon = item['next_taxon']
    closest = item['closest_genomes'][0]

    with open('CLOSEST_DISTANCE', 'w') as f:
      f.write(fmt_dist(closest['distance']))

    # output-writing function to reduce redunancy
    def write_output(file, search_item, column, empty_value):
      with open(file, 'w') as f:
        if search_item is None:
          f.write(empty_value)
        elif search_item[column] is None:
          f.write(empty_value)
        else:
          if str(empty_value) == str(fmt_dist(0)):
            f.write(fmt_dist(search_item[column]))
          else:
            # remove candidate sub-speciation from taxon name
            if column == 'name':
              gambit_name = search_item[column]
              gambit_name = re.sub(r'_[A-Za-z]+', '', gambit_name)  # This line is added to remove _X where X is any letter
              f.write(gambit_name)
            else:
              f.write(search_item[column])

    # Predicted taxon    
    write_output('PREDICTED_TAXON', predicted, 'name', 'NA')
    write_output('PREDICTED_TAXON_RANK', predicted, 'rank', 'NA')
    write_output('PREDICTED_TAXON_THRESHOLD', predicted, 'distance_threshold', fmt_dist(0))

    # Next taxon
    write_output('NEXT_TAXON', next_taxon, 'name', 'NA')
    write_output('NEXT_TAXON_RANK', next_taxon, 'rank', 'NA')
    write_output('NEXT_TAXON_THRESHOLD', next_taxon, 'distance_threshold', fmt_dist(0))

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
   
    # set merlin tags
    # The purpose of the merlin_tag output is for use as a trigger for organism-specific or taxon-specific workflows
    # One primary & important example is running NCBI amrfinderplus with the appropriate 'amrfinder --organism <organism>' option

    merlin_tag_designations = {"Escherichia" : "Escherichia", "Shigella" : "Escherichia", "Shigella sonnei" : "Shigella sonnei",
        "Klebsiella" : "Klebsiella", "Klebsiella pneumoniae" : "Klebsiella pneumoniae", "Klebsiella oxytoca" : "Klebsiella oxytoca", 
        "Klebsiella aerogenes" : "Klebsiella aerogenes", "Listeria" : "Listeria", "Salmonella" : "Salmonella", "Vibrio" : "Vibrio",
        "Vibrio cholerae" : "Vibrio cholerae"
    }

    try:
      merlin_tag = predicted['name']
      # remove candidate sub-speciation from merlin_tag
      merlin_tag = re.sub(r'_[A-Za-z]+', '', merlin_tag)  # This line is added to remove _X where X is any letter
    except:
      merlin_tag = "NA"

    # see if there is a reduced tag available (use Escherichia for Shigella flexneri)
    reduced_name = [val for key,val in merlin_tag_designations.items() if key in merlin_tag]

    if len(reduced_name) > 0: # if a reduced tag was identified, check if it should be used
      if (reduced_name[0] in merlin_tag_designations.keys()) and (merlin_tag not in merlin_tag_designations.keys()):
        merlin_tag = reduced_name[0]

    print(merlin_tag)
    with open('MERLIN_TAG', 'w') as merlin:
      merlin.write(merlin_tag)

    EOF
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
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: "~{cpu}"
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    maxRetries: 3
    preemptible: 1
  }
}
