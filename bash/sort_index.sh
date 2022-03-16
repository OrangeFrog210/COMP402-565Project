#!/bin/bash
module load samtools

# Assume that the user is in the sample director
# Assume that the user is in the directory of the individual (e.g. HG00615 in this case)
# For each loci (represented by a directory)

cd $1
for d1 in */
do
    # Move to the directory
    cd $1
    echo $d1
    cd $d1

    # Sort and index mapped.sam
    samtools view -b mapped.sam > mapped.bam
    samtools sort mapped.bam > mapped.sorted.bam
    samtools index mapped.sorted.bam

done

cd $1