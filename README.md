# RNA-seq_data_processing

## Download programs

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

## 1. Checking quality and trimming reads: FastQC and Trimmomatic

  Reads should be checked via FastQC, then trimmed of adapters and low quality sequences, then re-checked via FastQC. 
  
  Trimmomatic website and options: http://www.usadellab.org/cms/?page=trimmomatic
  
  **NOTE: If there are overrepresented sequences, add to the adapter file and rerun.**

1. Run 01a_RNA-seq_processing.py to set up .sh and .sub file for FastQC and Trimmomatic to be submitted to CHTC. 

Code and options:

> python 01a_RNA_seq_processing.py 
>   
Required:
>   
>   -fastq_list  a tab-delimited file with list of fastq files to process, for PE R1 in first column, R2 in second column
>   
>   -layout      PE: paired end, SE: single end
>   
>   -dir         the path to your staging directory on chtc
>   
Optional:
>   
>   -SRA         SRA sample ID (if starting from SRA file); default='NA'
>   
>   -threads     how many threads to be used for multiprocessing on chtc; default=4
>   
>   -adapters    adapter sequences to be trimmed; default='TruSeq3-PE.fa'
>   
>   -seedMis     seed mismatches: specifies the maximum mismatch count which will still allow a full match to be performed; default=2
>   
>   -pClipThres  palindromeClipThreshold: specifies how accurate the match between the two adapter ligated reads must be for PE palindrome read alignment; default=30
>   
>   -sClipThres  simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read; default=10
>   
>   -LEADING     Cut bases off the start of a read, if below a threshold quality- Specifies the minimum quality required to keep a base; default=3
>   
>   -TRAILING    Cut bases off the end of a read, if below a threshold quality- Specifies the minimum quality required to keep a base; default=3
>   
>   -windowSize  Specifies the number of bases to average quality across, cutting once the average quality within the window falls below a threshold; default=4
>   
>   -requiredQuality   Specifies the average quality required; default=20
>   
>   -MINLEN      Specifies the minimum length of reads to be kept; default=36
>   
>   -HEADCROP    Specifies number of bases at beginning of read to be removed; default=0
>   
>   -mem         Amount of memory asked to run job; default=100GB
>   
>   -disk        Amount of disk space asked to run job; default=100GB

Example running base code:

> python ~/Desktop/Github/RNA-seq_data_processing/01a_RNA-seq_processing.py -fastq_list file_list_initial.txt -layout PE -dir /staging/bmoore22

Output should be fastq_list.csv, fastqc_trimm_loop.sub, fastqc_trimm.sh

  See Examples folder for examples of these files.

2. Transfer output files to chtc- your home directory- and run:

    Make sure .fastq.gz files are in your staging directory

    Make sure fastqc_v0.11.9.zip and Trimmomatic-0.39.zip are in your home directory

    Submit to chtc:

    > condor_submit fastqc_trimm_loop.sub

    **NOTE: If you are rerunning with new sequences in the adapter file, make sure to add your new adapter file to the Trimmomatic folder and rezip.**

## 2. Mapping and counting reads: HiSAT2 and HTSeq

Read mapping is performed via HiSAT2 and then the reads are sorted via SamTools. Reads are then counted for each gene via HTSeq.

HiSAT2 manual: http://daehwankimlab.github.io/hisat2/manual/

Description of the Hisat2 sam format: http://daehwankimlab.github.io/hisat2/manual/#sam-output

HTSeq manual: https://htseq.readthedocs.io/en/release_0.11.1/count.html 

1. Run 01b_RNA-seq_processing.py to set up .sh file, .sub file, and file_list.csv for HiSAT2 and HTSeq on chtc.

Code and options:

> python 01b_RNA-seq_processing.py

Options:

Required:

> -trimmed_list   tab-delimited file with list of trimmed fastq files to input; for PE R1 in first column, R2 in second column
> -genome_seq     genome sequence fasta file
> -gff            gff3 file
> -layout         PE: paired end, SE: single end
> -dir            the path to your staging directory on chtc
> -base           base name for genome index

Optional:

> -threads        how many threads to be used; default=4
> -mem            mem asked for job running; default=60GB
> -disk           disk space asked for job running; default=100GB

Example running base code:

> python ~/Desktop/Github/RNA-seq_data_processing/01b_RNA-seq_processing.py -trimmed_list initial_trim_filelist.txt -genome_seq Joinvillea_ascendens_subsp.gabra.faa.mod.fa -gff Jascendensv1.1.gene.gff3 -layout PE -dir /staging/bmoore22 -base Jvasc

Output should be hisat2-htseq_filelist.csv, hisat2-htseq.sub, hisat2-htseq.sh

  See Examples folder for examples of these files.
  
2. Transfer output files to chtc- your home directory- and run:

    Make sure .fastq.trim.gz files are in your staging directory

    Make sure hisat2-2.2.1-Linux_x86_64.zip, samtools-1.15.1.tar.bz2, and packages2.tar.gz are in your home directory

    Submit to chtc:

    > condor_submit hisat2-htseq.sub


## 3. Combine the read count files

1. Transfer HTSeqCount files to your computer from chtc
2. In directory with HTSeqCount files, run 03_combine_read_counts.py:

> python 03_combine_read_counts.py

Output:

Read_counts.txt

This is your input for calculating either TPM/FPKM or for differential expression

## 4. Calculate expression (TPM and FPKM)

**Note: if you are only interested in differential expression, skip to step 5.**

1. Calculate transcript length (based on total number of exons for gene, using gff file):

Use:

> transcript_length.R

2. Make tissue-replicate file:

    a tab-delimited file with the .sorted_sam filename in column 1, and the tissue type in column 2
    
    see tissue-rep_list.txt in Example_data folder
    
3. Run 04_TPM-FPKM_calling.R

    INPUTs:
    
    Read_counts.txt file
    
    .transcript_length.txt file
    
    tissue-rep_list.txt file
    
## 5. Calculate Differential Expression

1. Use 05_DEseq.R
2. Try with sample data in Example_data folder: counts_data.csv first
3. INPUTS:

  Read_counts.txt file
  
  treatments list: currently to edit in code where each "treatment" is specified based on the column order. Replicates of same sample should be specified as the same treatment.


