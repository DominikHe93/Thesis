#!/bin/bash
#ran on unicluster
i=$(realpath $@)
paths=($(echo $i))
if [[ ! ${#paths[*]} -gt 0 ]];
	then
		echo -e "\nERROR: you need specify input fasta(s)\n"
		exit 2
fi
for p in ${paths[*]}; do
	if [[ ! -e "$p" ]];
		then
			echo -e "\nERROR: path \"$p\" does not exist\n"
			exit 2
	fi
done
workdir=$PWD
currdir=$(basename $PWD)
if [[ "$TMPDIR" == "" ]]; then
	TMPDIR="."
fi
echo "mktemp -p ${TMPDIR} -d XXXXXXX_rnammer${currdir}"
tempdir=$(mktemp -p ${TMPDIR} -d XXXXXXX_rnammer${currdir})
echo $tempdir


echo "cd ${tempdir}"
cd ${tempdir}
echo "extract_larger.py -if $i -co 500 -o temp_scaffolds.fasta"
extract_larger.py -if $i -co 500 -o temp_scaffolds.fasta

echo "grep -c ">" temp_scaffolds.fasta"
grep -c ">" temp_scaffolds.fasta

echo "subdivide_multifas.py -if temp_scaffolds.fasta -o subdivide -nf 4"
subdivide_multifas.py -if temp_scaffolds.fasta -o subdivide -nf 4

for inf in subdivide_*.fasta; do
	echo "rnammer -S bac-m lsu -f 23S_rRNA_${inf} $inf"
	rnammer -S bac -m lsu -f 23S_rRNA_${inf} $inf &
done 
echo "waiting for rnammer to finish..."
wait
echo "finished rnammer"
mkdir ${HOME}/delme
#echo "cp 23S_rRNA* ${HOME}/delme/."
#cp 23S_rRNA* ${HOME}/delme/.
echo "cat 23S_rRNA_*.fasta > ${workdir}/23S_rRNA_${currdir}.fasta || exit"
cat 23S_rRNA_*.fasta > ${workdir}/23S_rRNA_${currdir}.fasta || exit
echo "cd ${workdir}"
cd ${workdir}
echo "rm -rf ${tempdir}"
rm -rf ${tempdir}
#rm ${HOME}/delme/23S_rRNA*.fasta
