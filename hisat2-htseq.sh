#!/bin/bash
 
#unpack hisat2
unzip hisat2-2.2.1-Linux_x86_64.zip 
 
#change permissions
chmod -R 755 hisat2-2.2.1/
 
# load GCC module
export PATH
. /etc/profile.d/modules.sh
module load GCC/8.3.0
 
# test program is working
hisat2-2.2.1/hisat2
 
# copy file from staging
cp /staging/bmoore22/$1 ./
cp /staging/bmoore22/$2 ./
cp /staging/bmoore22/Jvasc_genome_index.zip ./
 
# unzip
gunzip $1
gunzip $2
unzip Jvasc_genome_index.zip
 
# copy index to working directory
cp Jvasc_genome_index/Jvasc* ./
 
# run hisat2
hisat2-2.2.1/hisat2 -x Jvasc -1 $3 -2 $4 -S $5 -q --dta --new-summary -p 4 --no-discordant
 
#rm trim files
rm *fastq.paired.trim*
 
# unzip samtools and install
tar xjf samtools-1.15.1.tar.bz2
chmod -R 777 samtools-1.15.1/
cd samtools-1.15.1/
./configure --disable-lzma
make
make install
cd ../
 
# run samtools sort by name
samtools-1.15.1/samtools sort -O sam -n $5 -o $6
 
# remove unsorted sam and index
rm $5
rm Jvasc_genome_index.zip
 
## untar your Python installation. Make sure you are using the right version!
tar -xzf python37.tar.gz
## untar packages
tar -xzf packages2.tar.gz
 
## export python and packages to working directory
export PATH=$PWD/python/bin:$PATH
export PYTHONPATH=$PWD/packages2
export HOME=$PWD
 
## run script
python3 -m HTSeq.scripts.count --format=sam -m union -s yes -t gene -i ID --nonunique=none -n 4 $6 Jascendensv1.1.gene.gff3 > $7
 
## remove sorted sam
rm $6

