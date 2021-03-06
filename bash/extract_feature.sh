#!/bin/bash
module load samtools
module load bedtools

cd $1
rm -f $2
# Assume that the user is in the sample director
# Assume that the user is in the directory of the individual (e.g. HG00615 in this case)
# For each loci (represented by a directory)
for d1 in */
do
	# Move to the directory
    cd $1
    cd $d1
	echo $d1
    d2=${d1::-1}

	# Extract the length of the REF genome from bam's header
	lengthREF_pre=$(samtools view -h mapped.sorted.bam | egrep "@SQ" | egrep "_genome")
	lengthREF=${lengthREF_pre#*LN:}
    echo "lengthREF: "
    echo "$lengthREF"
	# Calculate the right breakpoint's position
	ref_bkpt_r=$(($lengthREF - 500))
    #ref_bkpt_r=$(expr ${lengthREF} - 500);
    echo "ref_bkpt_r: "
	echo "$ref_bkpt_r"


	# Useful in later steps when counting the number of reads that spans a breakpoint or region
	genome_name="${d2}_genome"
	#echo "$genome_name"
	alt_name="${d2}_alternative"
	#echo "$alt_name"

    # Count the # of grey reads (MAPQ > 0)
    ## >> TEST
    COUNT_grey_REF=$(samtools view -h -q 1 mapped.sorted.bam | awk '{if ($3~"_genome") {print $0}}' | wc -l)
    COUNT_grey_ALT=$(samtools view -h -q 1 mapped.sorted.bam | awk '{if ($3~"_alternative") {print $0}}' | wc -l)
    COUNT_grey=$(samtools view -q 1 mapped.sorted.bam | wc -l) ## fixed. -h removed.

    
    # REF: # of grey reads which are mapped in proper pair and is MAPQ > 0 for both reads
    COUNT_2_REF=$(samtools view -f 0x2 -h mapped.sorted.bam | samtools view -h -q 1 | awk '{if ($3~"_genome" && $7=="=") {print $0}}' | wc -l)
    # ALT: # of grey reads which are mapped in proper pair and is MAPQ > 0 for both reads
    COUNT_2_ALT=$(samtools view -f 0x2 -h mapped.sorted.bam | samtools view -h -q 1 | awk '{if ($3~"_alternative" && $7=="=") {print $0}}' | wc -l)
    # Total: # of grey reads which are mapped in proper pair and is MAPQ > 0 for both reads
    COUNT_2=$(samtools view -f 0x2 -h mapped.sorted.bam | samtools view -h -q 1| awk '{if ($7=="=") {print $0}}' | wc -l)

    # Count overlap with the three breakpoints
	# Count the number of reads that overlap with the left breakpoint on REF (500bp position)
    
	count_overlap_ref_bkpt_l=$(samtools view -b mapped.sorted.bam ${genome_name}:500-500 | samtools view -q 1 | wc -l)
    #echo "$count_overlap_ref_bkpt_l"
	# Count the number of reads that overlap with the right breakpoint on REF (-500bp from the length of REF genome)
	count_overlap_ref_bkpt_r=$(samtools view -b mapped.sorted.bam ${genome_name}:${ref_bkpt_r}-${ref_bkpt_r} | samtools view -q 1 | wc -l)
    #echo "$count_overlap_ref_bkpt_r"
	# Count the number of reads that overlap with the breakpoint on ALT (500bp position)
	count_overlap_alt_bkpt=$(samtools view -b mapped.sorted.bam ${alt_name}:500-500 | samtools view -q 1 | wc -l)
    #echo "$count_overlap_alt_bkpt"
    
    # Extract header
    samtools view -H mapped.sorted.bam > header.txt
    samtools view -f 0x2 mapped.sorted.bam | awk '{if ($7=="=") {print $0}}' > body
    #samtools view -f 0x2 mapped.sorted.bam | awk '{if ($3~"_genome" && $7=="=") {print $0}}' > ref_body
    #samtools view -f 0x2 mapped.sorted.bam | awk '{if ($3~"_alternative" && $7=="=") {print $0}}' > alt_body
    
    cat header.txt body | samtools view -b - | bedtools bamtobed -i stdin | tr / $'\t' > temp.bed
    cut -f 4 temp.bed | sort | uniq -c | awk '$1==2 {print $2}' > 2reads.txt

    cat header.txt body | samtools view -b - | awk 'getline second {print $0"\t"second}' | awk '$1 == $8 {if ($6 > 0 && $13 > 0) {print $1"\t"$2"\t"$10"\t"$4"-"$11"\t"(($6+$13)/2)"\t\."}}' | awk '{if ($3<$2) {print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6} else {print $0}}' > reads.bed
    grep -f 2reads.txt temp.bed | sort -k4,4 -k2,2n | awk 'getline second {print $0"\t"second}' | awk '$1 == $8 {if ($6 > 0 && $13 > 0) {print $1"\t"$2"\t"$10"\t"$4"-"$11"\t"(($6+$13)/2)"\t\."}}' | awk '{if ($3<$2) {print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6} else {print $0}}' | sort -k2,2n > reads.bed
    #cat header.txt ref_body | samtools view -b - | bedtools bamtobed -i stdin | tr / $'\t' | sort -k4,4 | awk 'getline second {print $0"\t"second}' | awk '$1 == $8 {if ($6 > 0 && $13 > 0) {print $1"\t"$2"\t"$10"\t"$4"-"$12"\t"(($6+$13)/2)"\t\."}}' > ref_reads.bed
    #cat header.txt alt_body | samtools view -b - | bedtools bamtobed -i stdin | tr / $'\t' | sort -k4,4 | awk 'getline second {print $0"\t"second}' | awk '$1 == $8 {if ($6 > 0 && $13 > 0) {print $1"\t"$2"\t"$10"\t"$4"-"$12"\t"(($6+$13)/2)"\t\."}}' > alt_reads.bed
    
    # cat reads.bed
    # Preparing minigenome.bed
    rm -f minigenome.bed
    touch minigenome.bed
    echo -e "${genome_name}\t500\t501" > minigenome.bed
    ref_bkpt_r_1=$(($ref_bkpt_r + 1))
    echo "$ref_bkpt_r_1"
    echo -e "${genome_name}\t${ref_bkpt_r}\t${ref_bkpt_r_1}" >> minigenome.bed
    echo -e "${alt_name}\t500\t501" >> minigenome.bed
    head -3 minigenome.bed  # Check whether the content of minigenome.bed is correct
    # To get the counts
    # -wa: write the original entry in A for each overlap
    bedtools intersect -wa -b reads.bed -a minigenome.bed -loj | awk  '$5!="-1"' | awk '{print $1"_"$2"_"$3}'  | sort -k1,1 | uniq -c > count.txt

    join -a1 -11 -22 <(awk '{print $1"_"$2"_"$3}' minigenome.bed | sort -k1,1) <(sort -k2,2 count.txt) |  sed 's/ /\t/g' | awk '{if (NF==1) {print $0"\t0"} else {print $0} }' > count2.txt
    # <(): Subshell <()
    # a1: left outer join

    #bedtools intersect -wa -b reads.bed -a minigenome.bed -loj | cut -f 1-3 | sort | uniq -c > count.txt

    cat count2.txt
    
    alt_bkpr_frag=$(cat count2.txt | awk '{print $2}' | sed -n 1p)
    ref_bkpr_l_frag=$(cat count2.txt | awk '{print $2}' | sed -n 2p)
    ref_bkpr_r_frag=$(cat count2.txt | awk '{print $2}' | sed -n 3p)

	# Output features
	OUTPUT="$d2,$COUNT_grey_REF,$COUNT_grey_ALT,$COUNT_grey,$COUNT_2_REF,$COUNT_2_ALT,$COUNT_2,$count_overlap_ref_bkpt_l,$count_overlap_ref_bkpt_r,$count_overlap_alt_bkpt,$ref_bkpr_l_frag,$ref_bkpr_r_frag,$alt_bkpr_frag"
    cd $1
    echo $OUTPUT >> $2
    echo $OUTPUT
    echo "---------"
done
cd $1
head -10 $2