version 1.0

task tb_gene_coverage {
  input {
    File bamfile
    File bamindex
    String samplename
    Int min_depth = 10
    Int disk_size = 100
    String docker = "quay.io/staphb/samtools:1.15"
  }
  command <<<
    chr=$(samtools idxstats ~{bamfile} | cut -f 1 | head -1)

    # samtools outputs 3 columns; column 3 is the depth of coverage per nucleotide position, piped to awk to count the positions
    #  above min_depth, then wc -l counts them all
    gyrB=$(samtools depth -J -r "${chr}:5040-7467" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    gyrA=$(samtools depth -J -r "${chr}:7102-10018" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fgd1=$(samtools depth -J -r "${chr}:490583-491993" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    mshA=$(samtools depth -J -r "${chr}:575148-576990" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ccsA=$(samtools depth -J -r "${chr}:619691-621065" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rpoB=$(samtools depth -J -r "${chr}:759607-763525" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rpoC=$(samtools depth -J -r "${chr}:763170-767520" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    mmpL5=$(samtools depth -J -r "${chr}:775386-778680" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    mmpS5=$(samtools depth -J -r "${chr}:778277-779105" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    mmpR5=$(samtools depth -J -r "${chr}:778790-779687" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rpsL=$(samtools depth -J -r "${chr}:781360-782134" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rplC=$(samtools depth -J -r "${chr}:800609-801662" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fbiC=$(samtools depth -J -r "${chr}:1302731-1305701" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    Rv1258c=$(samtools depth -J -r "${chr}:1405881-1407540" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    embR=$(samtools depth -J -r "${chr}:1415981-1417547" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    atpE=$(samtools depth -J -r "${chr}:1460845-1461490" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rrs=$(samtools depth -J -r "${chr}:1471646-1473582" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rrl=$(samtools depth -J -r "${chr}:1473458-1476995" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fabG1=$(samtools depth -J -r "${chr}:1673148-1674383" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    inhA=$(samtools depth -J -r "${chr}:1673848-1675211" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rpsA=$(samtools depth -J -r "${chr}:1833342-1835187" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    tlyA=$(samtools depth -J -r "${chr}:1917740-1918946" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ndh=$(samtools depth -J -r "${chr}:2101451-2103242" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    katG=$(samtools depth -J -r "${chr}:2153689-2156570" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    PPE35=$(samtools depth -J -r "${chr}:2167449-2170812" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    Rv1979c=$(samtools depth -J -r "${chr}:2221519-2223364" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    pncA=$(samtools depth -J -r "${chr}:2288481-2290323" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    kasA=$(samtools depth -J -r "${chr}:2517915-2519565" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    eis=$(samtools depth -J -r "${chr}:2713924-2715586" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ahpC=$(samtools depth -J -r "${chr}:2725912-2726980" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    folC=$(samtools depth -J -r "${chr}:2745935-2747798" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    pepQ=$(samtools depth -J -r "${chr}:2859100-2860618" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ribD=$(samtools depth -J -r "${chr}:2986639-2987815" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    Rv2752c=$(samtools depth -J -r "${chr}:3064315-3066391" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    thyX=$(samtools depth -J -r "${chr}:3066993-3068161" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    thyA=$(samtools depth -J -r "${chr}:3073480-3074671" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ald=$(samtools depth -J -r "${chr}:3086620-3088135" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fbiD=$(samtools depth -J -r "${chr}:3338918-3339962" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    Rv3083=$(samtools depth -J -r "${chr}:3448304-3450191" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fprA=$(samtools depth -J -r "${chr}:3473807-3475577" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    whiB7=$(samtools depth -J -r "${chr}:3568201-3568879" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    Rv3236c=$(samtools depth -J -r "${chr}:3611759-3613316" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fbiA=$(samtools depth -J -r "${chr}:3640343-3641738" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    fbiB=$(samtools depth -J -r "${chr}:3641335-3643081" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    alr=$(samtools depth -J -r "${chr}:3839994-3841620" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    rpoA=$(samtools depth -J -r "${chr}:3877264-3878707" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ddn=$(samtools depth -J -r "${chr}:3986644-3987499" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    clpC1=$(samtools depth -J -r "${chr}:4037958-4040904" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    panD=$(samtools depth -J -r "${chr}:4043662-4044481" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    embC=$(samtools depth -J -r "${chr}:4239663-4243347" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    embA=$(samtools depth -J -r "${chr}:4243004-4246717" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    embB=$(samtools depth -J -r "${chr}:4246314-4250010" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    aftB=$(samtools depth -J -r "${chr}:4266753-4269036" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ubiA=$(samtools depth -J -r "${chr}:4268725-4270033" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ethA=$(samtools depth -J -r "${chr}:4325804-4330174" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    ethR=$(samtools depth -J -r "${chr}:4327349-4328399" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    whiB6=$(samtools depth -J -r "${chr}:4337971-4338721" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )
    gid=$(samtools depth -J -r "${chr}:4407328-4408476" ~{bamfile} | awk -F "\t" '{if ($3 > ~{min_depth}) print;}' | wc -l )

    # add one to gene lenth to compensate for subtraction
    gyrB_pc=$(python3 -c "print ( round ( ($gyrB / 2428 ) * 100, 2 ) )") 
    gyrA_pc=$(python3 -c "print ( round( ($gyrA / 2917 ) * 100, 2 ) )") 
    fgd1_pc=$(python3 -c "print ( round( ($fgd1 / 1411 ) * 100, 2 ) )")
    mshA_pc=$(python3 -c "print ( round( ($mshA / 1843 ) * 100, 2 ) )")
    ccsA_pc=$(python3 -c "print ( round( ($ccsA / 1375 ) * 100, 2 ) )")
    rpoB_pc=$(python3 -c "print ( round( ($rpoB / 3919 ) * 100, 2 ) )")
    rpoC_pc=$(python3 -c "print ( round( ($rpoC / 4351 ) * 100, 2 ) )")
    mmpL5_pc=$(python3 -c "print ( round( ($mmpL5 / 3295 ) * 100, 2 ) )")
    mmpS5_pc=$(python3 -c "print ( round( ($mmpS5 / 829 ) * 100, 2 ) )")
    mmpR5_pc=$(python3 -c "print ( round( ($mmpR5 / 898 ) * 100, 2 ) )")
    rpsL_pc=$(python3 -c "print ( round( ($rpsL / 775 ) * 100, 2 ) )")
    rplC_pc=$(python3 -c "print ( round( ($rplC / 1054 ) * 100, 2 ) )")
    fbiC_pc=$(python3 -c "print ( round( ($fbiC / 2971 ) * 100, 2 ) )")
    Rv1258c_pc=$(python3 -c "print ( round( ($Rv1258c / 1660 ) * 100, 2 ) )")
    embR_pc=$(python3 -c "print ( round( ($embR / 1567 ) * 100, 2 ) )")
    atpE_pc=$(python3 -c "print ( round( ($atpE / 646 ) * 100, 2 ) )")
    rrs_pc=$(python3 -c "print ( round( ($rrs / 1937 ) * 100, 2 ) )")
    rrl_pc=$(python3 -c "print ( round( ($rrl / 3538 ) * 100, 2 ) )")
    fabG1_pc=$(python3 -c "print ( round( ($fabG1 / 1236 ) * 100, 2 ) )")
    inhA_pc=$(python3 -c "print ( round( ($inhA / 1364 ) * 100, 2 ) )")
    rpsA_pc=$(python3 -c "print ( round( ($rpsA / 1846 ) * 100, 2 ) )")
    tlyA_pc=$(python3 -c "print ( round( ($tlyA / 1207 ) * 100, 2 ) )")
    ndh_pc=$(python3 -c "print ( round( ($ndh / 1792 ) * 100, 2 ) )")
    katG_pc=$(python3 -c "print ( round( ($katG / 2882 ) * 100, 2 ) )")
    PPE35_pc=$(python3 -c "print ( round( ($PPE35 / 3364 ) * 100, 2 ) )")
    Rv1979c_pc=$(python3 -c "print ( round( ($Rv1979c / 1846 ) * 100, 2 ) )")
    pncA_pc=$(python3 -c "print ( round( ($pncA / 1843 ) * 100, 2 ) )")
    kasA_pc=$(python3 -c "print ( round( ($kasA / 1651 ) * 100, 2 ) )")
    eis_pc=$(python3 -c "print ( round( ($eis / 1663 ) * 100, 2 ) )")
    ahpC_pc=$(python3 -c "print ( round( ($ahpC / 1069 ) * 100, 2 ) )")
    folC_pc=$(python3 -c "print ( round( ($folC / 1864 ) * 100, 2 ) )")
    pepQ_pc=$(python3 -c "print ( round( ($pepQ / 1519 ) * 100, 2 ) )")
    ribD_pc=$(python3 -c "print ( round( ($ribD / 1177 ) * 100, 2 ) )")
    Rv2752c_pc=$(python3 -c "print ( round( ($Rv2752c / 2077 ) * 100, 2 ) )")
    thyX_pc=$(python3 -c "print ( round( ($thyX / 1169 ) * 100, 2 ) )")
    thyA_pc=$(python3 -c "print ( round( ($thyA / 1192 ) * 100, 2 ) )")
    ald_pc=$(python3 -c "print ( round( ($ald / 1516 ) * 100, 2 ) )")
    fbiD_pc=$(python3 -c "print ( round( ($fbiD / 1045 ) * 100, 2 ) )")
    Rv3083_pc=$(python3 -c "print ( round( ($Rv3083 / 1888 ) * 100, 2 ) )")
    fprA_pc=$(python3 -c "print ( round( ($fprA / 1771 ) * 100, 2 ) )")
    whiB7_pc=$(python3 -c "print ( round( ($whiB7 / 679 ) * 100, 2 ) )")
    Rv3236c_pc=$(python3 -c "print ( round( ($Rv3236c / 1558 ) * 100, 2 ) )")
    fbiA_pc=$(python3 -c "print ( round( ($fbiA / 1396 ) * 100, 2 ) )")
    fbiB_pc=$(python3 -c "print ( round( ($fbiB / 1747 ) * 100, 2 ) )")
    alr_pc=$(python3 -c "print ( round( ($alr / 1627 ) * 100, 2 ) )")
    rpoA_pc=$(python3 -c "print ( round( ($rpoA / 1444 ) * 100, 2 ) )")
    ddn_pc=$(python3 -c "print ( round( ($ddn / 856 ) * 100, 2 ) )")
    clpC1_pc=$(python3 -c "print ( round( ($clpC1 / 2947 ) * 100, 2 ) )")
    panD_pc=$(python3 -c "print ( round( ($panD / 820 ) * 100, 2 ) )")
    embC_pc=$(python3 -c "print ( round( ($embC / 3685 ) * 100, 2 ) )")
    embA_pc=$(python3 -c "print ( round( ($embA / 3714 ) * 100, 2 ) )")
    embB_pc=$(python3 -c "print ( round( ($embB / 3697 ) * 100, 2 ) )")
    aftB_pc=$(python3 -c "print ( round( ($aftB / 2284 ) * 100, 2 ) )")
    ubiA_pc=$(python3 -c "print ( round( ($ubiA / 1309 ) * 100, 2 ) )")
    ethA_pc=$(python3 -c "print ( round( ($ethA / 4371 ) * 100, 2 ) )")
    ethR_pc=$(python3 -c "print ( round( ($ethR / 1051 ) * 100, 2 ) )")
    whiB6_pc=$(python3 -c "print ( round( ($whiB6 / 751 ) * 100, 2 ) )")
    gid_pc=$(python3 -c "print ( round( ($gid / 1149 ) * 100, 2 ) )")

    echo -e "#NOTE: THE VALUES BELOW ASSUME TBPROFILER (H37Rv) REFERENCE GENOME" > ~{samplename}.percent_gene_coverage.tsv
    echo -e "Gene\tPercent_Coverage" >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "gyrB\t"$gyrB_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "gyrA\t"$gyrA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fgd1\t"$fgd1_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "mshA\t"$mshA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ccsA\t"$ccsA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rpoB\t"$rpoB_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rpoC\t"$rpoC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "mmpL5\t"$mmpL5_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "mmpS5\t"$mmpS5_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv0678\t"$mmpR5_pc >> ~{samplename}.percent_gene_coverage.tsv # Rv0678 is the same as mmpR5
    echo -e "mmpR5\t"$mmpR5_pc >> ~{samplename}.percent_gene_coverage.tsv # Rv0678 is the same as mmpR5
    echo -e "rpsL\t"$rpsL_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rplC\t"$rplC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fbiC\t"$fbiC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv1258c\t"$Rv1258c_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "embR\t"$embR_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "atpE\t"$atpE_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rrs\t"$rrs_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rrl\t"$rrl_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fabG1\t"$fabG1_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "inhA\t"$inhA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rpsA\t"$rpsA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "tlyA\t"$tlyA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ndh\t"$ndh_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "katG\t"$katG_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "PPE35\t"$PPE35_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv1979c\t"$Rv1979c_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "pncA\t"$pncA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "kasA\t"$kasA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "eis\t"$eis_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ahpC\t"$ahpC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "folC\t"$folC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "pepQ\t"$pepQ_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ribD\t"$ribD_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv2752c\t"$Rv2752c_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "thyX\t"$thyX_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "thyA\t"$thyA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ald\t"$ald_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fbiD\t"$fbiD_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv3083\t"$Rv3083_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fprA\t"$fprA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "whiB7\t"$whiB7_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "Rv3236c\t"$Rv3236c_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fbiA\t"$fbiA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "fbiB\t"$fbiB_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "alr\t"$alr_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "rpoA\t"$rpoA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ddn\t"$ddn_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "clpC1\t"$clpC1_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "panD\t"$panD_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "embC\t"$embC_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "embA\t"$embA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "embB\t"$embB_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "aftB\t"$aftB_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ubiA\t"$ubiA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ethA\t"$ethA_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "ethR\t"$ethR_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "whiB6\t"$whiB6_pc >> ~{samplename}.percent_gene_coverage.tsv
    echo -e "gid\t"$gid_pc >> ~{samplename}.percent_gene_coverage.tsv

  >>>
  output {
    File tb_resistance_genes_percent_coverage = "~{samplename}.percent_gene_coverage.tsv"
  }
  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 2
    disks:  "local-disk " + disk_size + " SSD"
    disk: disk_size + " GB" 
    preemptible: 0
    maxRetries: 3
  }
}