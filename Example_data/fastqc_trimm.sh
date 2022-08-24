#!/bin/bash

#unpack fastqc
unzip fastqc_v0.11.9.zip 

#change permissions
chmod 755 FastQC/fastqc

# test program is working
FastQC/fastqc

# copy file from staging
cp /staging/bmoore22/$1 ./
cp /staging/bmoore22/$2 ./
# unzip file
gunzip $1
gunzip $2
# run fastqc
FastQC/fastqc -f fastq $3
FastQC/fastqc -f fastq $4
#unpack trimmomatic
unzip Trimmomatic-0.39.zip

# copy adapter sequences to working directory
cp Trimmomatic-0.39/adapters/TruSeq3-PE.fa ./

# run trimmmomatic on sample
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 4 $3 $4 $5 $6 $7 $8 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# rerun fastqc on new file
FastQC/fastqc -f fastq $5
FastQC/fastqc -f fastq $7

# gzip trimmed file
gzip *trim

# copy to staging
cp *trim.gz /staging/bmoore22/

# remove big files
rm *fastq
rm *trim.gz

# done- QC files will be transferred to home directory
