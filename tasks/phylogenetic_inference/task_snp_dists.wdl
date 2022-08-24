version 1.0

task snp_dists {
  input {
    File alignment
    String cluster_name
  }
  command <<<
    # date and version control
    date | tee DATE
    snp-dists -v | tee VERSION

    # create snp-dists matrix file
    snp-dists ~{alignment} > snp-dists-matrix.tsv 

    # create oredered snp-dists molten file 
    snp-dists -m ~{alignment} | awk '{print $NF,$0}' | sort -n | cut -f2- -d' ' > snp-dists-molten-ordered.tsv 

    # create list of isolates in order of SNP-dists 
    cut -f2 snp-dists-molten-ordered.tsv | awk '!seen[$0]++' > ordered-isolates.tsv

    python <<CODE
    # This script is based on Logan Fink's https://github.com/StaPH-B/CDPHE/blob/master/ordered_pairwise_matrix_generator_1.0.py
    # which takes in an ordered file and then returns a pairwise matrix based on that order.
    import sys, os
    cwd = os.getcwd()

    #Create a pairwise dictionary with all values possible, and create a list of the isolates
    pairwiseDict = {}
    isolates = []
    f = open("snp-dists-molten-ordered.tsv", 'r')
    flines  = f.readlines()
    for line in flines:
        print line
        line = line.strip('\n')
        linesplit = line.split('\t')
        pairwiseDict[linesplit[0], linesplit[1]] = linesplit[2]
        pairwiseDict[linesplit[1], linesplit[0]] = linesplit[2]
        if not linesplit[0] in isolates:
            isolates.append(linesplit[0])
        if not linesplit[1] in isolates:
            isolates.append(linesplit[1])
    f.close()

    #Take in the ordered list created by the user and create a list object
    orderedList = []
    g = open("ordered-isolates.tsv", 'r')
    glines = g.readlines()
    for line in glines:
        line = line.strip('\n')
        if str(line)[0] == "#":
            pass
        else:
            orderedList.append(line)
    g.close()

    # Seems inappropriate for snp-dists use case; removing for now
    #This is to make sure that if the names aren't exactly the same in the ordered file
    #(eg. doesn't include ".cleaned"), they will still work
    #newOrderedList = []
    #for item in orderedList:
    #    original = 1
    #    for isolate in isolates:
    #        if str(item) in str(isolate):
    #            newOrderedList.append(isolate)
    #            original = 0
    #        elif str(isolate) in str(item):
    #           newOrderedList.append(isolate)
    #           original = 0
    #        else:
    #            continue
    #    if original==1:
    #        newOrderedList.append(item)

    #Open the output file and write the first line
    z=open("~{cluster_name}_snp_distance_matrix.tsv", 'w')
    z.write(".")
    for item in orderedList:
        z.write('\t')
        z.write(str(item))
    z.write('\n')

    #Write the matrix using the ordered pairs dictionary
    for item1 in orderedList:
        z.write(str(item1))
        for item2 in orderedList:
            if item1 == item2:
                z.write('\t')
                z.write("-")
            else:
                z.write('\t')
                z.write(str(pairwiseDict[item1, item2]))
        z.write('\n')
    z.close()
    print "Matrix has been created in current directory as '~{cluster_name}_snp_distance_matrix.tsv.'"

    CODE
  >>>
  output {
    String date = read_string("DATE")
    String version = read_string("VERSION")
    File snp_matrix = "${cluster_name}_snp_distance_matrix.tsv"
    File snp_dists_molten_ordered = "snp-dists-molten-ordered.tsv"
  }
  runtime {
    docker: "quay.io/staphb/snp-dists:0.8.2"
    memory: "2 GB"
    cpu: 2
    disks: "local-disk 100 SSD"
    preemptible: 0
  }
}
