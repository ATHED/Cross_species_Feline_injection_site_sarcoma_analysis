#!/bin/bash

for name in lane2-s031-indexRPI1-ATCACG-Ela-C
do
   file_prefix=/nfs1/BIOMED/Ramsey_Lab/weiq/phase2_fastq_trimmed/${name}
   /local/cluster/STAR/bin/Linux_x86_64/STAR  --genomeDir /nfs1/BIOMED/Ramsey_Lab/weiq/RNAseq/RNAseqdata/reference_completegenome \
                                              --sjdbGTFfile /nfs1/BIOMED/Ramsey_Lab/weiq/RNAseq/RNAseqdata/reference_completegenome/Felis_catus.Felis_catus_6.2.87.gtf \
                                              --readFilesIn ${file_prefix}_R1_val_1.fq \
					                    ${file_prefix}_R2_val_2.fq \
                                              --outSAMtype BAM Unsorted SortedByCoordinate \
                                              --outFilterMultimapNmax 1 \
                                              --outSAMunmapped Within \
                                              --quantMode TranscriptomeSAM \
                                              --twopassMode Basic \
                                              --runThreadN 12 \
                                              --outFileNamePrefix  ${name}_star
done

