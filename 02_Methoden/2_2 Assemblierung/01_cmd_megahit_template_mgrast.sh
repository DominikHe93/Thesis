#!/bin/bash

sub_category="groundwater"
main_category="mg_rast"
prefix="$main_category""_""$sub_category"
thread=8
mem_flag="0.5"
min_kmer=21
max_kmer=99

min_contig=750
forw_reads=(Metagenome_*_pass_1.fastq.gz)
rev_reads=(Metagenome_*_pass_2.fastq.gz)
unpaired=(Metagenome_*_unpaired_*.fastq.gz *_pass.fastq.gz)
forw_arg=""
rev_arg=""
unpair_arg=""
if [[ -f "$forw_reads" || -f "${forw_reads[1]}" ]]
	then
		sep=" "
		index=0
		forw_arg="-1"
		rev_arg="-2"
		while [[ "$index" -lt "${#forw_reads[*]}" ]]
			do
				if [[ -f "${forw_reads[$index]}" ]]
				then
					forw_arg="${forw_arg}${sep}${forw_reads[$index]}"
					rev_arg="${rev_arg}${sep}${rev_reads[$index]}"
					sep=","
				fi
				index=$(($index+1))
		done
fi
if [[ -f "$unpaired" || -f "${unpaired[1]}" ]]
		then
			sep=" "
			index=0
			unpair_arg="-r"
			while [[ "$index" -lt "${#unpaired[*]}" ]]
			do
				if [[ -f "${unpaired[$index]}" ]]
				then
					unpair_arg="${unpair_arg}${sep}${unpaired[$index]}"
					sep=","
				fi
				index=$(($index+1))
			done
			
fi
echo "megahit --out-prefix $prefix -t $thread -m $mem_flag --k-min $min_kmer --k-max $max_kmer --k-step 10 $forw_arg $rev_arg $unpair_arg"
/usr/bin/time -v megahit --out-prefix $prefix -t $thread -m $mem_flag --k-min $min_kmer --k-max $max_kmer --kmin-1pass --min-contig-len $min_contig --k-step 10 $forw_arg $rev_arg $unpair_arg
echo "finished"
