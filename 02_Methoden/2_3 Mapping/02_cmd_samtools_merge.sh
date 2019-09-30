#!/bin/bash

kerne=4
for bamdir in mapping_*/; do
	sample=KITcoassembly_$(basename $bamdir| sed 's#mapping_##')
	#echo $sample
		# samtools merge -@ ANZAHL_KERNE -O FORMAT OUTPUT_DATEI INPUT_BAM_DATEIEN...
		echo "samtools merge -@ $kerne -O BAM ${bamdir}/merged_${sample}.bam ${bamdir}/${sample}*.bam && rm ${bamdir}/${sample}*.bam*"
		samtools merge -@ $kerne -O BAM ${bamdir}/merged_${sample}.bam ${bamdir}/${sample}*.bam && rm ${bamdir}/${sample}*.bam* #alle probenzugehörigen_bam dateien im ordner mergen und originale danach sofort löschen
	
		echo "samtools sort -@ $kerne -O BAM -o ${bamdir}/sorted_merged_${sample}.bam ${bamdir}/merged_${sample}.bam && rm ${bamdir}/merged_${sample}.bam"
		samtools sort -@ $kerne -O BAM -o ${bamdir}/sorted_merged_${sample}.bam ${bamdir}/merged_${sample}.bam && rm ${bamdir}/merged_${sample}.bam #eben germergte dateien sortieren. danach unsortierte originale sofort löschen

		echo "samtools index ${bamdir}/sorted_merged_${sample}.bam" 
		samtools index ${bamdir}/sorted_merged_${sample}.bam #sortierte gemergte ergebnisse inizieren

done


