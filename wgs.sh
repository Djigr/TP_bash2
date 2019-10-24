#!/bin/bash

mkdir $2
mkdir $2/index
mkdir $2/aligned

#Build index for database
#soft/bowtie2-build databases/all_genome.fasta databases/all_genome.fasta

#Build index for reads?
soft/bowtie2-build $1/*.fastq.gz $2/all_genome.fasta

#Alignment end-to-end
soft/bowtie2 -x databases/all_genome.fasta --end-to-end -1 $1/EchF_R1.fastq.gz -2 $1/EchF_R2.fastq.gz --very-fast -S $2/aligned/aligned.sam


#Convert from SAM to BAM
samtools view -h -b ../../$2/aligned/aligned.sam -o ../../$2/aligned/aligned.bam

#Sort the BAM
samtools sort $2/aligned/aligned.bam -o $2/aligned/alignsorted.bam

#Index the BAM
#samtools index $2/aligned/alignsorted.bam

#Extracting the counting
samtools idxstats $2/aligned/alignsorted.bam > $2/aligned/count.tsv

grep ">" databases/all_genome.fasta|cut -f 2 -d ">" >association.tsv

mkdir $2/megahit
#Question 3
soft/./megahit -1 $1/*1.fastq.gz -2 $1/*2.fastq.gz --k-list 21 --mem-flag 0 -o $2/megahiti

#Question 4 : Predict genes on contigs
soft/./prodigal -i $2/megahiti/final.contigs.fa -d $2/predict.fasta

#Question 5 : Select complete genes
sed "s:>:*\n>:g" $2/predict.fasta | sed -n "/partial=00/,/*/p"|grep -v "*" > $2/genes_full.fna

#Question 6 : blastn
soft/./blastn -query $2/genes_full.fna -evalue 1E-3 -db databases/resfinder.fna -out final.tsv -perc_identity 0.8 -qcov_hsp_perc 0.8 -outfmt "6 qseqid sseqid pident qcovs evalue bitscore" -best_hit_score_edge 1E-4

#Il y a visiblement des gènes de résistance parmi nos génomes. 

