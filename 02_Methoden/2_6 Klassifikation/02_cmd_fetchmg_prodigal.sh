#!/bin/bash

echo "fetchmg.pl"

echo "/data/ibg5_dominik/tools/fetchMG/fetchMG.pl -p -t 4 -m extraction /data/ibg5_dominik/project/results/prodigal/merged_assembly_totalprots.faa -o fetchmg_results"

/data/ibg5_dominik/tools/fetchMG/fetchMG.pl -p -t 4 -m extraction /data/ibg5_dominik/project/results/prodigal/merged_assembly_totalprots.faa -o fetchmg_results

echo "done"

#-p only proteinseq; -o create directory with results
