#!/bin/bash

#declare variables
database=/data/ibg5_dominik/project/results/merged_ncbi_mgrast_fasta/final_complete_merging_NCBIvsMGRAST.fasta

commonprefix="coassembly"

coupled1="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519753.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519753.3_pass_2.fastq.gz"

name1="mgm4519753.3"

input_string1="-c $coupled1 -p ${commonprefix}_${name1} -o mapping_$name1 -k"


coupled2="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519755.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519755.3_pass_2.fastq.gz"

name2="mgm4519755.3"

input_string2="-c $coupled2 -p ${commonprefix}_${name2} -o mapping_$name2 -k -K"


coupled3="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519757.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519757.3_pass_2.fastq.gz"

name3="mgm4519757.3"

input_string3="-c $coupled3 -p ${commonprefix}_${name3} -o mapping_$name3 -k -K"


coupled4="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519759.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519759.3_pass_2.fastq.gz"

name4="mgm4519759.3"

input_string4="-c $coupled4 -p ${commonprefix}_${name4} -o mapping_$name4 -k -K"


coupled5="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519761.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519761.3_pass_2.fastq.gz"

name5="mgm4519761.3"

input_string5="-c $coupled5 -p ${commonprefix}_${name5} -o mapping_$name5 -k -K"


coupled6="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519763.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519763.3_pass_2.fastq.gz"

name6="mgm4519763.3"

input_string6="-c $coupled6 -p ${commonprefix}_${name6} -o mapping_$name6 -k -K"


coupled7="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519765.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519765.3_pass_2.fastq.gz"

name7="mgm4519765.3"

input_string7="-c $coupled7 -p ${commonprefix}_${name7} -o mapping_$name7 -k -K"


coupled8="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519767.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519767.3_pass_2.fastq.gz"

name8="mgm4519767.3"

input_string8="-c $coupled8 -p ${commonprefix}_${name8} -o mapping_$name8 -k -K"


coupled9="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519769.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519769.3_pass_2.fastq.gz"

name9="mgm4519769.3"

input_string9="-c $coupled9 -p ${commonprefix}_${name9} -o mapping_$name9 -k -K"


coupled10="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519771.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519771.3_pass_2.fastq.gz"

name10="mgm4519771.3"

input_string10="-c $coupled10 -p ${commonprefix}_${name10} -o mapping_$name10 -k -K"


coupled11="/data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519773.3_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/mg_rast/190704_metagenome_fastq/Metagenome_mgm4519773.3_pass_2.fastq.gz"

name11="mgm4519773.3"

input_string11="-c $coupled11 -p ${commonprefix}_${name11} -o mapping_$name11 -k -K"


coupled12="/data/ibg5_dominik/project/results/downloaded_datasets/ncbi/190704/SRR6484324/SRR6484324_pass_1.fastq.gz /data/ibg5_dominik/project/results/downloaded_datasets/ncbi/190704/SRR6484324/SRR6484324_pass_2.fastq.gz"

name12="SRR6484324"

input_string12="-c $coupled12 -p ${commonprefix}_${name12} -o mapping_$name12 -k -K"


#actually do something


for i in "$input_string1" "$input_string2" "$input_string3" "$input_string4" "$input_string5" "$input_string6" "$input_string7" "$input_string8" "$input_string9" "$input_string10" "$input_string11" "$input_string12"; do
	echo "********************************************"
	
	echo "bamm make -t 4 -d $database $i"

	bamm make -t 4 -d $database $i
done

echo "done"



