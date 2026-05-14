#!/bin/bash

# Submit one PBS alignment job for each chromosome listed in chromssome_list.txt.
# Each submitted job receives the chromosome name through CHROM
# and runs the pipeline with set TOOL and ARM.

while IFS= read -r chrom; do
    # Skip empty lines in the chromosome list.
    [[ -z "$chrom" ]] && continue

    qsub -v CHROM="${chrom}",TOOL="winnowmap",ARM="p" alignment_job.pbs
done < chromssome_list.txt
