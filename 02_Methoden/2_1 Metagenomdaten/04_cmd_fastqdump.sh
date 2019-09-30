#!/bin/bash
if [ $# -ne 1 ];
	then
		echo -e "\n You need to supply EXACTLY one SRA-accession as argument"
		exit 2
fi
. activate py27
set -x
sraid=$1
fastq-dump --outdir ${sraid} --gzip --skip-technical --split-3 --read-filter pass --clip --defline-seq '@$sn[_$ac][_$si][$sg]/$ri' ${sraid}
