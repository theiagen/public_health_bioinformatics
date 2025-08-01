version 1.0

task cauris_cladetyper {
  input {
    File assembly_fasta
    String samplename
    Int kmer_size = 11
    Float max_distance = 0.1
    
    Int cpu = 8
    Int disk_size = 100
    String docker = "us-docker.pkg.dev/general-theiagen/staphb/gambit:1.0.0"
    Int memory = 16

    File ref_clade1 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1_GCA_002759435.3_Cand_auris_B8441_V3_genomic.fasta"
    String ref_clade1_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade1_GCA_002759435.3_Cand_auris_B8441_V3_genomic.gbff"
    File ref_clade2 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.fasta"
    String ref_clade2_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade2_GCA_003013715.2_ASM301371v2_genomic.gbff"
    File ref_clade3 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.fasta"
    String ref_clade3_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade3_GCF_002775015.1_Cand_auris_B11221_V1_genomic.gbff"
    File ref_clade4 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.fasta"
    String ref_clade4_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade4_GCA_003014415.1_Cand_auris_B11243_genomic.gbff"
    File ref_clade5 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.fasta"
    String ref_clade5_annotated = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade5_GCA_016809505.1_ASM1680950v1_genomic.gbff"
    File ref_clade6 = "gs://theiagen-public-resources-rp/reference_data/eukaryotic/candidozyma/Cauris_Clade6_GCA_032714025.1_ASM3271402v1_genomic.fasta"
    String? ref_clade6_annotated
    }
  command <<<
    set -euo pipefail

    gambit --version | tee VERSION

    # create gambit signature file for six clades + input assembly
    gambit signatures create -o my-signatures.h5 -k ~{kmer_size} -p ATGAC ~{ref_clade1} ~{ref_clade2} ~{ref_clade3} ~{ref_clade4} ~{ref_clade5} ~{ref_clade6} ~{assembly_fasta}
    # calculate distance matrix for all seven signatures
    gambit dist --qs my-signatures.h5 --square -o ~{samplename}_matrix.csv

    # parse matrix to see closest clade to input assembly
    ## sort by 8th column (distance against input sequence)
    ## take top two columns: header, header, top hit
    ## take bottom of these rows (top hit)
    ## grab only file name (top_clade)
    sort -k8 -t ',' "~{samplename}_matrix.csv" | head -2 | tail -n-1 | awk -F',' '{print$1}' > TOP_CLADE
    sort -k8 -t ',' "~{samplename}_matrix.csv" | head -2 | tail -n-1 | awk -F',' '{print$8}' > MIN_DISTANCE

    # create empty files for clade reference and clade type
    echo "None" > CLADEREF
    echo "" > CLADETYPE

    python3 <<CODE
    import os
    import sys
    # strip extensions in accord with gambit's approach and acquire the best matching reference
    def strip_extensions(filename, extensions):
        for ext in extensions:
            if filename.endswith(ext):
                return filename[:-len(ext)]
        return filename

    # determine if minimum observed gambit distance meets threshold
    with open("MIN_DISTANCE", 'r') as raw:
        try:
            min_distance = float(raw.read().strip())
        except ValueError:
            min_distance = 1
    if min_distance > ~{max_distance}:
        print("Input assembly does not match any clade")
        print(f"Minimum distance {min_distance} is greater than maximum allowed distance {~{max_distance}}")
        sys.exit(0)

    # identify the top clade in accord with gambit's naming scheme
    with open("TOP_CLADE", 'r') as raw:
        top_clade = raw.read().strip()
    gzip_exts = ('.gz')
    fa_exts = ('.fasta', '.fna', '.ffn', '.faa', '.frn', '.fa')
    ref2annotation = {"~{ref_clade1}": ["~{if defined(ref_clade1_annotated) then ref_clade1_annotated else 'None'}", "Clade1"], 
                      "~{ref_clade2}": ["~{if defined(ref_clade2_annotated) then ref_clade2_annotated else 'None'}", "Clade2"],
                      "~{ref_clade3}": ["~{if defined(ref_clade3_annotated) then ref_clade3_annotated else 'None'}", "Clade3"],
                      "~{ref_clade4}": ["~{if defined(ref_clade4_annotated) then ref_clade4_annotated else 'None'}", "Clade4"], 
                      "~{ref_clade5}": ["~{if defined(ref_clade5_annotated) then ref_clade5_annotated else 'None'}", "Clade5"], 
                      "~{ref_clade6}": ["~{if defined(ref_clade6_annotated) then ref_clade6_annotated else 'None'}", "Clade6"]}
    claderef, cladetype = "None", ""
    for ref, annotation in ref2annotation.items():
        if ref.endswith('/'):
            ref = ref[:-1]  # remove trailing slash if present
        edit_ref0 = os.path.basename(ref)
        edit_ref1 = strip_extensions(edit_ref0, gzip_exts)
        edit_ref2 = strip_extensions(edit_ref1, fa_exts)
        if top_clade == edit_ref2:
            claderef = annotation[0]
            cladetype = annotation[1]
            break
    
    # report top clade
    with open("CLADEREF", 'w') as claderef_file:
        claderef_file.write(claderef)
    with open("CLADETYPE", 'w') as cladetype_file:
        cladetype_file.write(cladetype)
    CODE
  >>>
  output {
    String gambit_version = read_string("VERSION")
    String gambit_cladetype = read_string("CLADETYPE")
    String annotated_reference = read_string("CLADEREF")
    String gambit_cladetyper_docker_image = docker
  }
  runtime {
    docker: docker
    memory: "~{memory} GB"
    cpu: cpu    
    disks: "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB"
    preemptible: 0
  }
}
