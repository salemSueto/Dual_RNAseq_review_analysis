#!/bin/bash

#The script takes a list of raw single-end fastq files. It performs pre-processing analysis 
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

#General
ngs_format=".fastq.gz" #RNA seq reads format
hostRef="rootDirectory/../reference/hsapiens" #Host reference directory
pathogenRef="rootDirectory/../reference/spneumoniae_D39" #Pathogen reference directory

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
##################


###Pre-processing#####
#Raw Data Quality Control
fastqc 0_raw_data/*$ngs_format -o 0_raw_data/qc
multiqc 0_raw_data/qc -o 0_raw_data/qc

#Raw Data Trimming
cd 0_raw_data
for f in *$ngs_format
do
    $trimmomatic SE $f ../1_trimmed_data/${f%%$ngs_format}"_trim.$ngs_format" ILLUMINACLIP:$adapter:$illumina_clip SLIDINGWINDOW:$sw MINLEN:$minlen;
done
cd ../

#Trimmed Data Quality Control
$fastqc 1_trimmed_data/*$ngs_extension -o 1_trimmed_data/qc
$multiqc 1_trimmed_data/qc -o 1_trimmed_data/qc


#####STAR mapping#####
#Host
cd 1_trimmed_data
for f in *$ngs_extension
do 
    $star --genomeDir $hostRef/star --readFilesIn $f --quantMode TranscriptomeSAM GeneCounts --outFilterMismatchNoverLmax 0.05 --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate    --outFileNamePrefix ../2_star_host/${f%%$ngs_extension}"."; 
done
cd ../
mv 2_star_host/*ReadsPerGene* 2_star_host/htseq

#Pathogen
cd 1_trimmed_data

for f in *$ngs_extension
do 
    $star --genomeDir $pathogenRef/star --readFilesIn $f --readFilesCommand zcat --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.05 --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix ../3_star_pathogen/${f%%$ngs_extension}"."; 
done
cd ../
mv 3_star_pathogen/*ReadsPerGene* 3_star_pathogen/htseq


#####RSEM quantification#####
#Host
cd  2_star_host
for f in *toTranscriptome*
do 
    $rsem-calculate-expression --bam --no-bam-output --estimate-rspd --forward-prob $fp $f $hostRef/rsem/index rsem/${f%%.toTranscriptome.out.bam};
done
cd ..

#Pathogen
cd 3_star_pathogen
for f in *toTranscriptome*
do 
    $rsem-calculate-expression --bam --no-bam-output --estimate-rspd --forward-prob $fp $f $pathogenRef/rsem/index rsem/${f%%.toTranscriptome.out.bam};
done
cd ..


#####Salmon#####
#Host
cd 1_trimmed_data
for f *$ngs_extension
do
    $salmon quant -i $hostRef/salmon -l A -r $f -p 8 --validateMappings -o ../4_salmon_host/${f%%$ngs_extension} --numBootstraps 10;
done
cd ..

#Pathogen
cd 1_trimmed_data
for f in *$ngs_extension
do
    $salmon quant -i $pathogenRef/salmon -l A -r $f -p 8 --validateMappings -o ../5_salmon_pathogen/${f%%$ngs_extension} --numBootstraps 10;
done
cd ../


#####Kallisto#####
#Host
cd 1_trimmed_data
for f in *$ngs_extension 
do
    $kallisto quant -i $hostRef/kallisto/index -l $avReadLength -o ../6_kallisto_host/${f%%$ngs_extension} -b 10 --single -s $sdReadLen $f;
done
cd ../

#Pathogen
cd 1_trimmed_data
for f in *$ngs_extension
do
    $kallisto quant -i $pathogenRef/kallisto/index -l $avReadLength -o ../7_kallisto_pathogen/${f%%$ngs_extension} -b 10 --single  -s $sdReadLen $f;
done
cd ..
