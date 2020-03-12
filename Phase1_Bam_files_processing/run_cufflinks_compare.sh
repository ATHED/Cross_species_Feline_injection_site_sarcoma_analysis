#!/bin/bash

# save a log
exec &> run_logfile_cuffcompare.mod.txt

/local/cluster/cufflinks-2.2.1.Linux_x86_64/cuffcompare -V -r /nfs1/BIOMED/Ramsey_Lab/weiq/RNAseq/RNAseqdata/reference_completegenome/Felis_catus.Felis_catus_6.2.87.mod.gtf transcript.mod.gtf -o fiss_cufflinks_mod


