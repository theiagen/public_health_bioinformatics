version 1.0

task aa_subs {
  meta {
    description: "Calls amino acid mutations from a given alignment file"
  }
  input {
    File alignment
    String protein_name
    String? subtype
    # place holder for organism variable, might be used in the future for organisms other than flu
    # String? organism
    String docker = "us-docker.pkg.dev/general-theiagen/broadinstitute/viral-core:2.1.33"
    Int disk_size = 100
    Int cpu = 2
  }
  command <<<
  if [[ "~{protein_name}" == "NA" && "~{subtype}" == "H1N1" ]] ; then
    reference_id="CY121682.1"
  elif [[ "~{protein_name}" == "NA" && "~{subtype}" == "H3N2" ]] ; then
    reference_id="CY114383.1"
  elif [[ "~{protein_name}" == "HA" && "~{subtype}" == "H1N1" ]] ; then
    reference_id="MW626062.1"
  elif [[ "~{protein_name}" == "HA" && "~{subtype}" == "H3N2" ]] ; then
    reference_id="EPI_ISL_1563628"
  elif [[ "~{protein_name}" == "PA" ]] ; then
    reference_id="CY121685.1"
  elif [[ "~{protein_name}" == "PB1" ]] ; then
    reference_id="CY121686.1"
  elif [[ "~{protein_name}" == "PB2" ]] ; then
    reference_id="CY121687.1"
  elif [[ "~{protein_name}" == "M2" && "~{subtype}" == "H1N1" ]] ; then
    reference_id="NC_026431.1"
  elif [[ "~{protein_name}" == "M2" && "~{subtype}" == "H3N2" ]] ; then
    reference_id="NC_007367.1"
  else
    echo "Please enter a valid protein name such as HA, NA, PA, PB1, PB2, or M2 and a valid subtype such as H1N1 or H3N2."
  fi

  python3 <<CODE
  from Bio import AlignIO
  from Bio.Seq import Seq
  
  # Load the alignment file
  #alignment_file_path = 'sequences_aln.fasta'  # Replace with your actual file path
  alignment = AlignIO.read("~{alignment}", "fasta")
  
  # Identify the reference sequence
  reference_id = "${reference_id}"
  reference_sequence_record = next((record for record in alignment if reference_id in record.id), None)
  
  if reference_sequence_record is None:
    raise ValueError(f"Reference sequence with ID ${reference_id} not found in the alignment.")
  
  # Find the position of the first and last non-gap character in the reference sequence
  first_base_pos = next((pos for pos, base in enumerate(reference_sequence_record.seq) if base != '-'), None)
  last_base_pos = next((pos for pos, base in enumerate(reversed(reference_sequence_record.seq)) if base != '-'), None)
  
  # If the alignment is to be trimmed at the end, calculate the correct end position
  if last_base_pos is not None:
    last_base_pos = len(reference_sequence_record.seq) - last_base_pos
  
  # Ensure that the reference sequence is properly defined
  if first_base_pos is None or last_base_pos is None:
    raise ValueError("Could not find valid start and end positions in the reference sequence.")
  
  # Trim the reference sequence to remove leading and trailing gap columns
  reference_sequence = reference_sequence_record.seq[first_base_pos:last_base_pos]
  
  # Function to replace ambiguous nucleotides with 'N'
  def replace_ambiguous_nucleotides(sequence):
    return sequence.upper().replace('-', 'N')
  
  # Function to translate the sequence and handle 'NNN' codons
  def safe_translate(sequence):
    return str(Seq(replace_ambiguous_nucleotides(sequence)).translate())
  
  # Function to find amino acid mutations with protein prefix
  def find_amino_acid_mutations_with_protein_prefix(ref_protein, query_protein, protein_prefix):
    mutations = []
    for (position, (ref_aa, query_aa)) in enumerate(zip(ref_protein, query_protein), start=1):
      if ref_aa != query_aa and query_aa != '*':  # Exclude stop codons (*)
        mutations.append(f"{protein_prefix}:{ref_aa}{position}{query_aa}")
    return mutations
  
  # Translate the reference sequence to amino acids
  reference_protein = safe_translate(reference_sequence)
  
  # Open a TSV file for writing the results
  output_file_path = 'aa_changes.tsv' 
  with open(output_file_path, 'w') as output_file:
    # Write the header to the TSV file
    output_file.write("samplename\taa_changes\n")

    # Iterate over each sequence in the alignment
    for record in alignment:
      if record.id != reference_sequence_record.id:  # Exclude the reference sequence
        # Trim the query sequence to remove leading and trailing gap columns
        query_sequence = record.seq[first_base_pos:last_base_pos]
        # Translate the query sequence to amino acids
        query_protein = safe_translate(query_sequence)
        # Find mutations
        mutations = find_amino_acid_mutations_with_protein_prefix(reference_protein, query_protein, "~{protein_name}")
        mutations_string = ','.join(mutations) if mutations else 'None'
        if mutations:  # Only print if there are mutations
          print(f"{record.id} mutations: {' '.join(mutations)}")
          output_file.write(f"{record.id}\t{mutations_string}\n")
        else:
          print(f"{record.id} has no mutations relative to the reference.")
          output_file.write(f"{record.id}\tno mutations relative to the reference.\n")
  CODE
  >>>
  output {
    File aa_changes_tsv = "aa_changes.tsv"
  }
  runtime {
    docker: "~{docker}"
    memory: "8 GB"
    cpu: cpu
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}

task antiviral_mutations_parser {
  meta {
    description: "Parsing the output files from aa_subs task to identify mutations associated with influenza antiviral resistsnce."
  }
  input {
    File mutations_tsv
    String docker = "python:slim"
    Int disk_size = 50
    # User defined aa substitutions
    String? antiviral_aa_subs
  }
  command <<<
    touch TAMIFLU_AASUBS
    touch A_315675_AASUBS
    touch AMANTADINE_AASUBS
    touch COMPOUND_367_AASUBS
    touch FAVIPIRAVIR_AASUBS
    touch FLUDASE_AASUBS
    touch L_742_001_AASUBS
    touch LANINAMIVIR_AASUBS
    touch PERAMIVIR_AASUBS
    touch PIMODIVIR_AASUBS
    touch RIMANTADINE_AASUBS
    touch OSELTAMIVIR_AASUBS
    touch XOFLUZA_AASUBS
    touch ZANAMIVIR_AASUBS
    
    # Parse outputs using python3
    python3 <<CODE
    import csv
    import codecs

    # list of aa substitutions linked with resistance to various antivirals
    A_315675_aa_subs = ["NA:A272V", "NA:D151E", "NA:D151N", "NA:E119D", "NA:E119G", "NA:E119V", "NA:H274Y", "NA:H275Y", "NA:N146K", "NA:Q136K", "NA:R150K", "NA:R152K", "NA:R289K", "NA:R292K", "NA:S219T", "NA:V116A"]
    A_315675_aa_subs = A_315675_aa_subs + "~{antiviral_aa_subs}".split(',')

    amantadine_aa_subs = ["M2:L26F", "M2:L26I", "M2:V27A", "M2:A30S", "M2:A30T", "M2:S31N", "M2:G34E"]
    amantadine_aa_subs = amantadine_aa_subs + "~{antiviral_aa_subs}".split(',')

    compound_367_aa_subs = ["PB1:H456P"]
    compound_367_aa_subs = compound_367_aa_subs + "~{antiviral_aa_subs}".split(',')

    favipiravir_aa_subs = ["NA:H275Y", "NA:H274Y"]
    favipiravir_aa_subs = favipiravir_aa_subs + "~{antiviral_aa_subs}".split(',')

    fludase_aa_subs = ["HA:G125R", "HA:S125T", "HA:G126R", "HA:S136T", "HA:G137R", "HA:S186I", "NA:W438L"]
    fludase_aa_subs = fludase_aa_subs + "~{antiviral_aa_subs}".split(',')

    L_742_001_aa_subs = ["PA:T20A", "PA:H41A", "PA:G81F", "PA:G81T", "PA:G81V", "PA:F105S", "PA:I120T"]
    L_742_001_aa_subs = L_742_001_aa_subs + "~{antiviral_aa_subs}".split(',')

    laninamivir_aa_subs = ["HA:K130E", "HA:V132A", "HA:K133E", "HA:K134E", "HA:V135A", "HA:S138A", "HA:K153E", "HA:G155E", "HA:K156E", "HA:G158E", "HA:P194L", "HA:D222G", "HA:D225G", "NA:G104E", "NA:E105K", "NA:G109E", "NA:E110K", "NA:E117G", "NA:E117A", "NA:E119G", "NA:E119A", "NA:E119K", "NA:E119V", "NA:E119D", "NA:Q136R", "NA:Q136K", "NA:P139S", "NA:G140R", "NA:P141S", "NA:N142S", "NA:G142R", "NA:G147E", "NA:T148I", "NA:R150K", "NA:D151N", "NA:R152K", "NA:I219K", "NA:I219R", "NA:I222K", "NA:I222R", "NA:I223K", "NA:I223R", "NA:S246R", "NA:S247R", "NA:H274Y", "NA:H275Y", "NA:R289K", "NA:R292K", "NA:R293K", "NA:Q313R", "NA:I427T"]
    laninamivir_aa_subs = laninamivir_aa_subs + "~{antiviral_aa_subs}".split(',')

    peramivir_aa_subs = ["NA:A136S", "NA:A138S", "NA:A200T", "NA:A201T", "NA:A245T", "NA:A246T", "NA:A390E", "NA:A395E", "NA:D151E", "NA:D151G", "NA:D151N", "NA:D197E", "NA:D197N", "NA:D197Y", "NA:D198E", "NA:D198G", "NA:D198N", "NA:D198Y", "NA:D199G", "NA:D199N", "NA:D429G", "NA:D432G", "NA:E99A", "NA:E99D", "NA:E99G", "NA:E105E", "NA:E105K", "NA:E110E", "NA:E110K", "NA:E117A", "NA:E117D", "NA:E117G", "NA:E117V", "NA:E119A", "NA:E119D", "NA:E119G", "NA:E119I", "NA:E119V", "NA:G104E", "NA:G109E", "NA:G140R", "NA:G142R", "NA:G145E", "NA:G145R", "NA:G147E", "NA:G147R", "NA:H255Y", "NA:H271Y", "NA:H273Y", "NA:H274H", "NA:H274Y", "NA:H275H", "NA:H275Y", "NA:I117V", "NA:I120V", "NA:I122V", "NA:I203M", "NA:I203V", "NA:I219R", "NA:I219T", "NA:I221I", "NA:I221L", "NA:I221N", "NA:I221T", "NA:I221V", "NA:I222I", "NA:I222K", "NA:I222L", "NA:I222M", "NA:I222N", "NA:I222R", "NA:I222T", "NA:I222V", "NA:I223K", "NA:I223M", "NA:I223R", "NA:I223V", "NA:K358E", "NA:K359E", "NA:K360E", "NA:N142S", "NA:N144K", "NA:N146K", "NA:N275S", "NA:N294S", "NA:N295S", "NA:P139S", "NA:P141S", "NA:Q136K", "NA:Q136R", "NA:Q138K", "NA:Q138R", "NA:Q140K", "NA:Q140R", "NA:R150K", "NA:R152K", "NA:R221Q", "NA:R222Q", "NA:R289K", "NA:R290K", "NA:R292K", "NA:R293K", "NA:R371K", "NA:R374K", "NA:S246N", "NA:S246R", "NA:S247N", "NA:S247R", "NA:S334N", "NA:S336N", "NA:T106P", "NA:T110P", "NA:T111P", "NA:T146I", "NA:T148I", "NA:Y142H", "NA:Y144H", "NA:Y155H"]
    peramivir_aa_subs = peramivir_aa_subs + "~{antiviral_aa_subs}".split(',')

    pimodivir_aa_subs = ["PB2:Q306H", "PB2:S324I", "PB2:S324N", "PB2:S324R", "PB2:F404Y", "PB2:M431I", "PB2:N510T"]
    pimodivir_aa_subs = pimodivir_aa_subs + "~{antiviral_aa_subs}".split(',')
    
    rimantadine_aa_subs = ["M2:L26F", "M2:L26I", "M2:V27A", "M2:A30S", "M2:A30T", "M2:S31N", "M2:G34E"]
    rimantadine_aa_subs = rimantadine_aa_subs + "~{antiviral_aa_subs}".split(',')

    oseltamivir_aa_subs = ["HA:A28T", "HA:D222G", "HA:D225G", "HA:G155E", "HA:G158E", "HA:K130E", "HA:K133E", "HA:K134E", "HA:K140E", "HA:K144E", "HA:K153E", "HA:K156E", "HA:K234Q", "HA:K238Q", "HA:P194L", "HA:R188K", "HA:R192K", "HA:R453M", "HA:S138A", "HA:T82K", "HA:T92K", "HA:V132A", "HA:V135A", "NA:A200T", "NA:A201A", "NA:A201T", "NA:A245T", "NA:A246T", "NA:A272V", "NA:A390E", "NA:A395E", "NA:D151D", "NA:D151E", "NA:D151G", "NA:D151N", "NA:D179G", "NA:D179N", "NA:D197E", "NA:D197N", "NA:D197Y", "NA:D198D", "NA:D198E", "NA:D198G", "NA:D198N", "NA:D198Y", "NA:D199E", "NA:D199G", "NA:D199N", "NA:D213G", "NA:D214G", "NA:D344N", "NA:D345N", "NA:D432G", "NA:E59G", "NA:E99A", "NA:E99D", "NA:E99G", "NA:E115V", "NA:E117A", "NA:E117D", "NA:E117G", "NA:E117V", "NA:E118V", "NA:E119A", "NA:E119D", "NA:E119E", "NA:E119G", "NA:E119I", "NA:E119K", "NA:E119V", "NA:E222V", "NA:E258Q", "NA:E272Q", "NA:E273Q", "NA:E276D", "NA:G104E", "NA:G108E", "NA:G109E", "NA:G140R", "NA:G142R", "NA:G145R", "NA:G147E", "NA:G147R", "NA:G320E", "NA:H101L", "NA:H255Y", "NA:H271Y", "NA:H273Y", "NA:H274H", "NA:H274N", "NA:H274Y", "NA:H275H", "NA:H275Y", "NA:H276Y", "NA:H277Y", "NA:H439P", "NA:H439R", "NA:I97V", "NA:I117M", "NA:I117V", "NA:I203L", "NA:I203M", "NA:I203R", "NA:I203T", "NA:I203V", "NA:I219K", "NA:I219L", "NA:I219M", "NA:I219R", "NA:I219T", "NA:I221L", "NA:I221N", "NA:I221T", "NA:I222K", "NA:I222L", "NA:I222M", "NA:I222N", "NA:I222R", "NA:I222T", "NA:I222V", "NA:I223K", "NA:I223L", "NA:I223M", "NA:I223R", "NA:I223T", "NA:I223V", "NA:I294V", "NA:I314V", "NA:I427T", "NA:K130N", "NA:K150N", "NA:K273Q", "NA:M372K", "NA:M375K", "NA:N44S", "NA:N46S", "NA:N142S", "NA:N144K", "NA:N146K", "NA:N169S", "NA:N199S", "NA:N200S", "NA:N220K", "NA:N221K", "NA:N275S", "NA:N294S", "NA:N295S", "NA:N325K", "NA:N329K", "NA:N368K", "NA:N369K", "NA:N386K", "NA:N390K", "NA:P139S", "NA:P141S", "NA:Q116L", "NA:Q136K", "NA:Q136L", "NA:Q313R", "NA:R136K", "NA:R150K", "NA:R151W", "NA:R152K", "NA:R152W", "NA:R193G", "NA:R194G", "NA:R221Q", "NA:R222Q", "NA:R224K", "NA:R289K", "NA:R290K", "NA:R292K", "NA:R293K", "NA:R371K", "NA:R374K", "NA:S219T", "NA:S227N", "NA:S245N", "NA:S246G", "NA:S246N", "NA:S246R", "NA:S247G", "NA:S247N", "NA:S247P", "NA:S247R", "NA:S331R", "NA:S334N", "NA:S336N", "NA:T146K", "NA:T146P", "NA:T148I", "NA:T156I", "NA:T157I", "NA:V95A", "NA:V96A", "NA:V116A", "NA:V215I", "NA:V233M", "NA:V234M", "NA:V240I", "NA:V241I", "NA:Y142H", "NA:Y144H", "NA:Y155H"]
    oseltamivir_aa_subs = oseltamivir_aa_subs + "~{antiviral_aa_subs}".split(',')

    xofluza_aa_subs = ["PA:I38F", "PA:I38M", "PA:I38T"]
    xofluza_aa_subs = xofluza_aa_subs + "~{antiviral_aa_subs}".split(',')

    zanamivir_aa_subs = ["HA:A28T", "HA:K140E", "HA:K144E", "HA:K162E", "HA:K165E", "HA:L317P", "HA:L320P", "HA:R453M", "NA:A200T", "NA:A201A", "NA:A201T", "NA:A243T", "NA:A245T", "NA:A246T", "NA:D151A", "NA:D151D", "NA:D151E", "NA:D151G", "NA:D151N", "NA:D151V", "NA:D179G", "NA:D197E", "NA:D197N", "NA:D197Y", "NA:D198D", "NA:D198E", "NA:D198G", "NA:D198N", "NA:D198Y", "NA:D199N", "NA:E99A", "NA:E99D", "NA:E99G", "NA:E105K", "NA:E110K", "NA:E117A", "NA:E117D", "NA:E117G", "NA:E118D", "NA:E118G", "NA:E119A", "NA:E119D", "NA:E119G", "NA:E119I", "NA:E119V", "NA:E276D", "NA:G104E", "NA:G109E", "NA:G140R", "NA:G142R", "NA:G402S", "NA:G407S", "NA:H255Y", "NA:H274H", "NA:H274Y", "NA:H275Y", "NA:I117M", "NA:I117V", "NA:I219R", "NA:I221L", "NA:I221N", "NA:I222K", "NA:I222L", "NA:I222N", "NA:I222R", "NA:I222V", "NA:I223K", "NA:I223R", "NA:I223V", "NA:I275T", "NA:I427T", "NA:N142S", "NA:N144K", "NA:N146K", "NA:N275S", "NA:N294S", "NA:N295S", "NA:N329K", "NA:P124T", "NA:P136T", "NA:P139S", "NA:P141S", "NA:P165L", "NA:Q116L", "NA:Q136K", "NA:Q136L", "NA:Q136Q", "NA:Q136R", "NA:Q138R", "NA:Q140R", "NA:Q313R", "NA:R150K", "NA:R151W", "NA:R152K", "NA:R152W", "NA:R224K", "NA:R289K", "NA:R290K", "NA:R292K", "NA:R293K", "NA:R371K", "NA:R374K", "NA:S110F", "NA:S245N", "NA:S246R", "NA:S247P", "NA:S247R", "NA:S249G", "NA:S250G", "NA:T40A", "NA:T43A", "NA:T106I", "NA:T110I", "NA:T111I", "NA:T148I", "NA:T148K", "NA:V96A", "NA:V116A", "NA:Y155H"]
    zanamivir_aa_subs = zanamivir_aa_subs + "~{antiviral_aa_subs}".split(',')

    # returns intersection between aa substitutions in input.tsv and known antiviral associated aa substitutions
    def intersection(present_mutations, known_mutations):
      return list(set(present_mutations) & set(known_mutations))

    # writes the intersection of aa substitutions associated with each antiviral to a file
    def write_to_file(filename, aa_mutations, known_mutations):
      mutations = intersection(aa_mutations, known_mutations.split(','))
      with codecs.open(filename, 'wt') as file:
        file.write(", ".join(mutations))

    # read in aa substitutions from input.tsv file
    with codecs.open("~{mutations_tsv}",'r') as tsv_file:
      tsv_reader = csv.reader(tsv_file, delimiter="\t")
      mutation_data = list(tsv_reader)

      # prevent failure if no mutations are found
      if len(mutation_data) == 1:
        mutation_data.append(['NA']*len(mutation_data[0]))
      tsv_dict = dict(zip(mutation_data[0], mutation_data[1]))

      with codecs.open("ALL_AASUBS", 'wt') as All_AA_Subs:
        nc_aa_subs = tsv_dict["aa_changes"]
        if nc_aa_subs == "no mutations relative to the reference.":
          nc_aa_subs = "NA"
        else:
          write_to_file("A_315675_AASUBS", A_315675_aa_subs, nc_aa_subs)
          write_to_file("AMANTADINE_AASUBS", amantadine_aa_subs, nc_aa_subs)
          write_to_file("COMPOUND_367_AASUBS", compound_367_aa_subs, nc_aa_subs)
          write_to_file("FAVIPIRAVIR_AASUBS", favipiravir_aa_subs, nc_aa_subs)
          write_to_file("FLUDASE_AASUBS", fludase_aa_subs, nc_aa_subs)
          write_to_file("L_742_001_AASUBS", L_742_001_aa_subs, nc_aa_subs)
          write_to_file("LANINAMIVIR_AASUBS", laninamivir_aa_subs, nc_aa_subs)
          write_to_file("PERAMIVIR_AASUBS", peramivir_aa_subs, nc_aa_subs)
          write_to_file("PIMODIVIR_AASUBS", pimodivir_aa_subs, nc_aa_subs)
          write_to_file("OSELTAMIVIR_AASUBS", tamiflu_aa_subs, nc_aa_subs)
          write_to_file("RIMANTADINE_AASUBS", rimantadine_aa_subs, nc_aa_subs)
          write_to_file("XOFLUZA_AASUBS", xofluza_aa_subs, nc_aa_subs)
          write_to_file("ZANAMIVIR_AASUBS", zanamivir_aa_subs, nc_aa_subs)

        All_AA_Subs.write(nc_aa_subs)
    CODE
  >>>
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
  output {
    String A_315675_aa_subs = read_string("A_315675_AASUBS")
    String amantadine_aa_subs = read_string("AMANTADINE_AASUBS")
    String compound_367_aa_subs = read_string("COMPOUND_367_AASUBS")
    String favipiravir_aa_subs = read_string("FAVIPIRAVIR_AASUBS")
    String fludase_aa_subs = read_string("FLUDASE_AASUBS")
    String L_742_001_aa_subs = read_string("L_742_001_AASUBS")
    String laninamivir_aa_subs = read_string("LANINAMIVIR_AASUBS")
    String peramivir_aa_subs = read_string("PERAMIVIR_AASUBS")
    String pimodivir_aa_subs = read_string("PIMODIVIR_AASUBS")
    String oseltamivir_aa_subs = read_string("TAMIFLU_AASUBS")
    String rimantadine_aa_subs = read_string("RIMANTADINE_AASUBS")
    String xofluza_aa_subs = read_string("XOFLUZA_AASUBS")
    String zanamivir_aa_subs = read_string("ZANAMIVIR_AASUBS")
  }
}

task serialization {
  input{ 
    Array[String] flu_A_315675_resistance_array = [""]
    Array[String] flu_amantadine_resistance_array = [""]
    Array[String] flu_compound_367_resistance_array = [""]
    Array[String] flu_favipiravir_resistance_array = [""]
    Array[String] flu_fludase_resistance_array = [""]
    Array[String] flu_L_742_001_resistance_array = [""]
    Array[String] flu_laninamivir_resistance_array = [""]
    Array[String] flu_peramivir_resistance_array = [""]
    Array[String] flu_pimodivir_resistance_array = [""]
    Array[String] flu_oseltamivir_resistance_array = [""]
    Array[String] flu_rimantadine_resistance_array = [""]
    Array[String] flu_xofluza_resistance_array = [""]
    Array[String] flu_zanamivir_resistance_array = [""]

    Int disk_size = 10
    String docker = "python:slim"
  }
  command <<<
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_A_315675_resistance_array}'.split(',') if mutation != '']))" > flu_A_315675_resistance.txt

    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_amantadine_resistance_array}'.split(',') if mutation != '']))" > flu_amantadine_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_compound_367_resistance_array}'.split(',') if mutation != '']))" > flu_compound_367_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_favipiravir_resistance_array}'.split(',') if mutation != '']))" > flu_favipiravir_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_fludase_resistance_array}'.split(',') if mutation != '']))" > flu_fludase_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_L_742_001_resistance_array}'.split(',') if mutation != '']))" > flu_L_742_001_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_laninamivir_resistance_array}'.split(',') if mutation != '']))" > flu_laninamivir_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_peramivir_resistance_array}'.split(',') if mutation != '']))" > flu_peramivir_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_pimodivir_resistance_array}'.split(',') if mutation != '']))" > flu_pimodivir_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_oseltamivir_resistance_array}'.split(',') if mutation != '']))" > flu_oseltamivir_resistance.txt

    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_rimantadine_resistance_array}'.split(',') if mutation != '']))" > flu_rimantadine_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_xofluza_resistance_array}'.split(',') if mutation != '']))" > flu_xofluza_resistance.txt
    
    python3 -c "print(', '.join([mutation for mutation in '~{sep=',' flu_zanamivir_resistance_array}'.split(',') if mutation != '']))" > flu_zanamivir_resistance.txt
  >>>
  output{
    String flu_A_315675_resistance = read_string("flu_A_315675_resistance.txt")
    String flu_amantadine_resistance = read_string("flu_amantadine_resistance.txt")
    String flu_compound_367_resistance = read_string("flu_compound_367_resistance.txt")
    String flu_favipiravir_resistance = read_string("flu_favipiravir_resistance.txt")
    String flu_fludase_resistance = read_string("flu_fludase_resistance.txt")
    String flu_L_742_001_resistance = read_string("flu_L_742_001_resistance.txt")
    String flu_laninamivir_resistance = read_string("flu_laninamivir_resistance.txt")
    String flu_peramivir_resistance = read_string("flu_peramivir_resistance.txt")
    String flu_pimodivir_resistance = read_string("flu_pimodivir_resistance.txt")
    String flu_oseltamivir_resistance = read_string("flu_tamiflu_resistance.txt")
    String flu_rimantadine_resistance = read_string("flu_rimantadine_resistance.txt")
    String flu_xofluza_resistance = read_string("flu_xofluza_resistance.txt")
    String flu_zanamivir_resistance = read_string("flu_zanamivir_resistance.txt")
  }
  runtime {
    docker: "~{docker}"
    memory: "4 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    dx_instance_type: "mem1_ssd1_v2_x2"
    maxRetries: 3
  }
}