#!/bin/bash

#The script takes a list of raw paired-end fastq files. It performs pre-processing analysis 
#and mapping with STAR, RSEM, Salmon, and Kallisto.


#####User input#####

#Tool's path
fastqc="fastqc"
multiqc="multiqc"
trimmomatic="trimmomatic"
star="star"
rsem-calculate-expression="rsem-calculate-expression"
salmon="salmon"
kallisto="kallisto"

#General information
ngs_format=".fastq.gz" #RNA seq reads format
hostRef="rootDirectory/../reference/hsapiens"
pathogenRef="rootDirectory/../reference/spneumoniae_D39"

#Pre-processing: trimmomatic
adapter="rootDirectory/../TruSeq3-SE.fa" #Adapter fasta file absolute path
illumina_clip="2:30:10" #Cut adapter and other illumina-specific sequences from the read
sw="5:20" #Sliding window & Cut within the window once the average quality falls below the threshold.
minlen="50" #Drop the read if it is below a specified length

#RSEM options
fp="1.0" #RSEM's forward_prob setting
avReadLength="75"

#Kallisto
avReadLength="75" #Average Read Length
sdReadLen="20" #Estimated standard deviation of fragment length

#Sample
sampleList=("sample_01" "sample_02" "sample_03" "sample_04" "sample_05" "sample_06" "sample_07" "sample_08" "sample_09" "sample_10")
##################


#####Pre-processing#####
#Raw Data Quality Control
$fastqc 0_raw_data/*$ngs_format -o 0_raw_data/qc
$multiqc 0_raw_data/qc -o 0_raw_data/qc

#Raw Data Trimming
cd 0_raw_data
for s in "${sampleList[@]}"
do
    $trimmomatic PE -p 8 $s*1* $s*2* $s*1"_paired$ngs_format" $s*1"_unpaired$ngs_format" $s*2"_paired"$ngs_format $s*2"_unpaired"$ngs_format ILLUMINACLIP:$adapter:$illumina_clip SLIDINGWINDOW:$sliding_window MINLEN:$minlen;
done
cd ../
mv 0_raw_data/*paired* 1_trimmed_data

#Trimmed Data Quality Control
$fastqc 1_trimmed_data/*_paired* -o 1_trimmed_data/qc
$multiqc 1_trimmed_data/qc -o 1_trimmed_data/qc


#####STAR mapping#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    $star --genomeDir $hostRef --readFilesIn $s*1_paired* $s*2_paired --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ../2_star_host/$s".";
done
cd ../
mv 2_star_host/*ReadsPerGene* 2_star_host/htseq

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do 
    $star --genomeDir $pathogenRef --readFilesIn $s*1_paired* $s*2_paired* --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --alignIntronMax 1 --outFileNamePrefix ../3_star_pathogen/$s".";
done
cd ../
mv 3_star_pathogen/*ReadsPerGene* 3_star_pathogen/htseq


#####RSEM quantification#####
#Host
cd 2_star_host
for s in *toTranscriptome*
do
    s0=${s%%Aligned.toTranscriptome.out.bam}
    $rsem-calculate-expression -paired-end --forward-prob $fp -p 8 --bam --no-bam-output --estimate-rspd $s $hostRef/rsem/index rsem/$s0; 
done
cd ../

#Pathogen
cd 2_star_pathogen
for s in *toTranscriptome*
do
    s0=${s%%.Aligned.toTranscriptome.out.bam}
    $rsem-calculate-expression --paired-end --forward-prob $fp -p 8 --bam --no-bam-output --estimate-rspd $f $pathogenRef/rsem/index rsem/$s0; 
done
cd ../


#####Salmon#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    $salmon quant -i $hostRef/salmon -l A -1 $s*1_paired* -2 $s*2_paired* -p 8 --validateMappings -o ../4_salmon_host/$s --numBootstraps 10;
done
cd ../

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    $salmon quant -i $pathogenRef/salmon -l A -1 $s*1_paired* -2 $s*2_paired* -p 8 --validateMappings -o ../5_salmon_pathogen/$s --numBootstraps 10;
done
cd ../


#####Kallisto#####
#Host
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    kallisto quant -i $hostRef/index -o ../6_kallisto_host/$s -b 10 $s*1_paired* $s*2_paired*;
done
cd ../

#Pathogen
cd 1_trimmed_data
for s in "${sampleList[@]}"
do
    $kallisto quant -i $pathogenRef/index -o ../7_kallisto_pathogen/$s -b 10 $s*1_paired* $s*2_paired*;
done
cd ../
