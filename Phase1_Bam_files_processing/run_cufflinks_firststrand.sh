#!/bin/bash

# save a log
exec &> run_logfile_cufflink_mod_1ststrand.txt

/local/cluster/cufflinks-2.2.1.Linux_x86_64/cufflinks -p 12 --library-type fr-firststrand -g /nfs1/BIOMED/Ramsey_Lab/weiq/RNAseq/RNAseqdata/reference_completegenome/Felis_catus.Felis_catus_6.2.87.mod.gtf -o transcript.mod.1ststrand.gtf out.bam


