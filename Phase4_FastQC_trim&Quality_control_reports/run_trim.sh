#!/bin/bash

trimgalore=/nfs1/BIOMED/Ramsey_Lab/software/trim_galore_zip
for name in `ls ../phase2_fastq/*.fastq.gz | cut -f1-2 -d_ | uniq`
do
   ${trimgalore}/trim_galore --paired \
		--gzip \
		--fastqc_args "--outdir ../phase2_fastq_trimmed_fastqc/" \
		${name}_R1.fastq.gz \
		${name}_R2.fastq.gz 2>>run_trim.log
done
