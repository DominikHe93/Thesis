#!/bin/bash

prodigal -q -a merged_assembly_totalprots.faa -i ~/project/results/merged_ncbi_mgrast_fasta/final_complete_merging_NCBIvsMGRAST.fasta -o /dev/null -p meta
