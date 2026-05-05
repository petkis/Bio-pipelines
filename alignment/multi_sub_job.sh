#!/bin/bash

while IFS= read -r chrom; do
    qsub -v CHROM="${chrom}",TOOL="winnowmap" paternal_gen_job.pbs
done < chroms_list.txt
