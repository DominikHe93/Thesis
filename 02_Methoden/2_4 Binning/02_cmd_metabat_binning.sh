#!/bin/bash

assembly=/data/ibg5_dominik/project/results/merged_ncbi_mgrast_fasta/final_complete_merging_NCBIvsMGRAST.fasta

sample1=/data/ibg5_dominik/project/results/mapping/mapping_SRR6484324/coassembly_SRR6484324final_complete_merging_NCBIvsMGRAST.SRR6484324_pass_1.bam

sample2=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519773.3/coassembly_mgm4519773.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519773.3_pass_1.bam

sample3=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519771.3/coassembly_mgm4519771.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519771.3_pass_1.bam

sample4=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519769.3/coassembly_mgm4519769.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519769.3_pass_1.bam

sample5=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519767.3/coassembly_mgm4519767.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519767.3_pass_1.bam

sample6=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519765.3/coassembly_mgm4519765.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519765.3_pass_1.bam

sample7=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519763.3/coassembly_mgm4519763.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519763.3_pass_1.bam

sample8=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519761.3/coassembly_mgm4519761.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519761.3_pass_1.bam

sample9=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519759.3/coassembly_mgm4519759.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519759.3_pass_1.bam

sample10=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519757.3/coassembly_mgm4519757.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519757.3_pass_1.bam

sample11=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519755.3/coassembly_mgm4519755.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519755.3_pass_1.bam

sample12=/data/ibg5_dominik/project/results/mapping/mapping_mgm4519753.3/coassembly_mgm4519753.3final_complete_merging_NCBIvsMGRAST.Metagenome_mgm4519753.3_pass_1.bam

echo "runMetaBat.sh -t 4 -m 1500 ${assembly} ${sample1} ${sample2} ${sample3} ${sample4} ${sample5} ${sample6} ${sample6} ${sample7} ${sample8} ${sample9} ${sample10} ${sample11} ${sample12}"

runMetaBat.sh -t 4 -m 1500 ${assembly} ${sample1} ${sample2} ${sample3} ${sample4} ${sample5} ${sample6} ${sample6} ${sample7} ${sample8} ${sample9} ${sample10} ${sample11} ${sample12}

echo "done"
