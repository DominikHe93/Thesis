#!/bin/bash

contig=/data/ibg5_dominik/project/results/merged_ncbi_mgrast_fasta/final_complete_merging_NCBIvsMGRAST.fasta

abund1=/data/ibg5_dominik/project/results/maxbin/myabundfile1.tab
abund2=/data/ibg5_dominik/project/results/maxbin/myabundfile2.tab
abund3=/data/ibg5_dominik/project/results/maxbin/myabundfile3.tab
abund4=/data/ibg5_dominik/project/results/maxbin/myabundfile4.tab
abund5=/data/ibg5_dominik/project/results/maxbin/myabundfile5.tab
abund6=/data/ibg5_dominik/project/results/maxbin/myabundfile6.tab
abund7=/data/ibg5_dominik/project/results/maxbin/myabundfile7.tab
abund8=/data/ibg5_dominik/project/results/maxbin/myabundfile8.tab
abund9=/data/ibg5_dominik/project/results/maxbin/myabundfile9.tab
abund10=/data/ibg5_dominik/project/results/maxbin/myabundfile10.tab
abund11=/data/ibg5_dominik/project/results/maxbin/myabundfile11.tab
abund12=/data/ibg5_dominik/project/results/maxbin/myabundfile12.tab

output=maxbinned_coassembly.fa

echo "run_MaxBin.pl -thread 4 -abund1 ${abund1} -abund2 ${abund2} -abund3 ${abund3} -abund4 ${abund4} -abund5 ${abund5} -abund6 ${abund6} -abund7 ${abund7} -abund8 ${abund8} -abund9 ${abund9} -abund10 ${abund10} -abund11 ${abund11} -abund12 ${abund12} -plotmarker -contig ${contig} -out ${output}"


run_MaxBin.pl -thread 4 -abund1 ${abund1} -abund2 ${abund2} -abund3 ${abund3} -abund4 ${abund4} -abund5 ${abund5} -abund6 ${abund6} -abund7 ${abund7} -abund8 ${abund8} -abund9 ${abund9} -abund10 ${abund10} -abund11 ${abund11} -abund12 ${abund12} -plotmarker -contig ${contig} -out ${output}


echo "done"
