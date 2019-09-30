#!/bin/bash

echo "gtdbtk metabat"

echo "gtdbtk classify_wf --cpus 4 -x .fa --genome_dir /data/ibg5_dominik/project/results/metabat/final_complete_merging_NCBIvsMGRAST.fasta.metabat-bins1500/ --out_dir gtdbtk_metabat_output"

gtdbtk classify_wf --cpus 4 -x .fa --genome_dir /data/ibg5_dominik/project/results/metabat/final_complete_merging_NCBIvsMGRAST.fasta.metabat-bins1500/ --out_dir gtdbtk_metabat_output

echo "done"

echo "*************************************"

echo "gtdbtk maxbin"

gtdbtk classify_wf --cpus 4 -x .fasta --genome_dir /data/ibg5_dominik/project/results/maxbin/maxbinned_coassembly.fasta/  --out_dir gtdbtk_maxbin_output

echo "done"

echo "**************************************"
