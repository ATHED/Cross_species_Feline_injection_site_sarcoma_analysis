#!/bin/bash
gtf_file=/nfs1/BIOMED/Ramsey_Lab/weiq/RNAseq/RNAseqdata/reference_completegenome/Felis_catus.Felis_catus_6.2.87.gtf
bamfilenames=`ls /nfs1/BIOMED/Ramsey_Lab/weiq/phase1_bam/*sortedByCoord.out.bam`
fileprefix=`/bin/basename ${bamfilename} .bam`     
/nfs1/BIOMED/Ramsey_Lab/software/bin/featureCounts -Q 3 \
						  -a ${gtf_file} \
						  -o counts.txt \
						  -T 10 \
						  ${bamfilenames}  2>get_counts_fc.log

