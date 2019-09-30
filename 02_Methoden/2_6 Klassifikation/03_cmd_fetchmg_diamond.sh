#!/bin/bash
set -x


#echo "diamond blastx --query /data/ibg5_dominik/project/results/fetchmg/fetchmg_results/markers.fasta --db /data/db/diamond_nr_20190115.dmnd --threats 4 --outfmt 6"

diamond blastp --query /data/ibg5_dominik/project/results/fetchmg/fetchmg_results/markers.fasta --db /data/db/diamond_nr_20190115.dmnd --threads 4 --outfmt 6 --out fetchmg_blasted


# --outfmt 6 BLAST tabular ; --query input data; --db database ; --out outputfilname; outfmt outputfileformat



echo "done"
