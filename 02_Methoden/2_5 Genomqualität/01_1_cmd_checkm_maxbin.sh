#!/bin/bash

input=/data/ibg5_dominik/project/results/maxbin/maxbinned_coassembly.fasta/

output=./taxonomy_wf_bacteria/

rank=domain

taxon=Bacteria

echo "taxonomy_wf"

echo "checkm taxonomy_wf -t 4 -x .fasta -f checkm_taxonomy_wf_overview_qa.tab --tab_table ${rank} ${taxon} ${input} ${output}"

checkm taxonomy_wf -t 4 -x .fasta -f checkm_taxonomy_wf_overview_qa.tab --tab_table ${rank} ${taxon} ${input} ${output}

echo "taxonomy done"

echo "*****************************************"

echo "lineage_wf"

output2=./lineage_wf/

echo "checkm lineage_wf -t 4 -x .fasta -f checkm_lineage_wf_overview_qa.tab --tab_table ${input} ${output2}"

checkm lineage_wf -t 4 -x .fasta -f checkm_lineage_wf_overview_qa.tab --tab_table ${input} ${output2}

echo "lineage done"

echo "*****************************************"

echo "checkm tree_qa"

echo "checkm tree_qa -o 2 -f tree_qa.tab --tab_table ./lineage_wf/"

checkm tree_qa -o 2 -f tree_qa.tab --tab_table ./lineage_wf/

echo "checkm tree_qa done"

echo "finish"
