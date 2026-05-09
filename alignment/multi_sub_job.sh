#!/bin/bash

while IFS= read -r chrom; do
    qsub -v CHROM="${chrom}",TOOL="winnowmap" alignment_job.pbs
done < chromssome_list.txt
