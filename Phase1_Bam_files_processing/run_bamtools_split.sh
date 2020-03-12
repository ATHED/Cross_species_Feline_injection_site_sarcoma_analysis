#!/bin/bash

# save a log
exec &> run_logfile_bamtools_split.txt

/nfs1/BIOMED/Ramsey_Lab/weiq/phase1_2_bam_fusion/bamtools/bin/bamtools-2.4.1 split -in out.bam -reference


