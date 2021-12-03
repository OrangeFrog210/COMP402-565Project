#!/bin/bash

mkdir ~/work/ae/867b940ddb7a3627d4b853c9f0b292/genotyping/IGV

for d1 in */ ; do
	samtools view -b mapped.sam > mapped.bam
	samtools sort mapped.bam > mapped.sorted.bam
	samtools index mapped.sorted.bam
	mkdir ~/work/ae/867b940ddb7a3627d4b853c9f0b292/genotyping/"$d1"
	mv mapped.sorted.bam ~/work/ae/867b940ddb7a3627d4b853c9f0b292/genotyping/IGV/"$d1"/"$("$d1" + "mapped.sorted.bam")"
done

for d2 in ~/work/ae/867b940ddb7a3627d4b853c9f0b292/genotyping/locusAlleles/; do
	cd "$d2"
	cp alleles.fa ~/work/ae/867b940ddb7a3627d4b853c9f0b292/genotyping/IGV/"$d2"/"$("$d2" + "alleles.fa")" 
done