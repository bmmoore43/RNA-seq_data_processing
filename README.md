# RNA-seq_data_processing

## download programs

**NOTE: if running on CHTC, these programs are already in the maeda_share drive**

Need:

FastQC

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Trimmomatic

download binary release zip: http://www.usadellab.org/cms/?page=trimmomatic

HISAT2

http://daehwankimlab.github.io/hisat2/download/

SAMtools

http://www.htslib.org/download/

Python HTSeq

**NOTE: Python should already be installed**

on command line:
> pip install HTSeq

SRA-Toolkit (if starting from an SRA file)

https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
  #### Note that if the fastq-dump doesn't work, run the code below to save the config for sra toolkit
  > vdb-config --interactive  


## build the genome index
> hisat2-build Pvirgatum_516_v5.0.fa Pvirgatum_516_v5.0.fa_gi

## You can apply the code 01_RNA-seq_processing.py to do the jobs automatically. An example command line for this code is as follow, which will produce a .sh file to be submitted to the slurm queue

> python 01_RNA_seq_processing.py -SRA SRR14066076 -genome_seq Pvirgatum_516_v5.0.fa -gff Pvirgatum_516_v5.1.gene.gff3 -layout PE -workdir /mnt/scratch/peipeiw/Data_for_Kenia/For_RNA_seq_pipeline -trim y -adapters all_PE_adapters.fa -base Pvirgatum_516_v5

#### Note: one step is conducted to keep only uniquely mapped reads which had Mapping quality of HISAT2 = 60. You may want to change the script 02_keep_reads_with_quality_60_and_unique_mapping.py/02_keep_reads_with_quality_60_and_unique_mapping_SE.py if that is not what you want. Description of the Hisat2 sam format: http://daehwankimlab.github.io/hisat2/manual/#sam-output
#### Note: if you want to keep only uniquely mapped reads, please check whether the read IDs for fastq_1 and fastq_2 are the same for pairs in the paired end data. If not the same, your output would be problematic.

## combine the read count files
> python 03_combine_read_counts.py

## call the TPM, which will be done in R
> Rscript 04_TPM_calling.r Read_counts.txt Pvirgatum_516_v5.1_transcript_length.txt
