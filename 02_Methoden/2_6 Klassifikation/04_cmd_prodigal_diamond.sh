#!/bin/bash

echo "diamond blastx --query /data/ibg5_dominik/project/results/prodigal/merged_assembly_totalprots.faa --db /data/db/diamond_nr_20190115.dmnd --threats 4 --outfmt 6"

diamond blastp --query /data/ibg5_dominik/project/results/prodigal/merged_assembly_totalprots.faa --db /data/db/diamond_nr_20190115.dmnd --threads 4 --outfmt 6 --out prodigal_blasted


# --outfmt 6 BLAST tabular ; --query input data; --db database ; --out outputfilname; outfmt outputfileformat



echo "done"
