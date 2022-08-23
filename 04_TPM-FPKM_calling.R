if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scater")
BiocManager::install("preprocessCore")

library(scater)
library('stringr')
library(devtools)
library(Biobase)
library(preprocessCore)

args = commandArgs(TRUE)
RC_file <- args[1] # your read counts file
TL_file <- args[2] # transcript length file

RC_file <- "Read_counts.txt"
TL_file <- "Jascendensv1.1.gene.gff3.transcript_length.txt"

Exp <- read.table(RC_file,head=T,sep='\t',stringsAsFactors=F,row.names=1)
dat <- read.table(TL_file,head=T,sep='\t',stringsAsFactors=F,row.names=1)
#Expression <- cbind(Exp[order(row.names(Exp)),],'Length'=dat[order(row.names(dat)),])
#Expression <- as.data.frame(Expression)
Expression <- merge(Exp, dat,
                          by = 'row.names', all = TRUE)
#Expression <- Expression[6:29125,]
#samp2 <- Expression[,-1]
#rownames(samp2) <- Expression[,1]
Len <- Expression$transcript_length
tpm_table <- calculateTPM(as.matrix(Expression[,2:(ncol(Expression)-1)]), Len)
tpm_table <- cbind(gene=Expression[,1],tpm_table)
#rownames(tpm_table) <- Expression[,1]
#write.table(tpm_table,'HTSeqCount_9-month-old-internode-2_sam3_sorted_quality_60_unique.sam.out_TPM.txt',row.names=T,sep='\t',quote=F)
#colnames(tpm_table)
fpkm_table <- calculateFPKM(as.matrix(Expression[,2:(ncol(Expression)-1)]), Len)
#rownames(fpkm_table) <- Expression[,1]
fpkm_table <- cbind(gene=Expression[,1],fpkm_table)
#write.table(fpkm_table,'HTSeqCount_9-month-old-internode-2_sam3_sorted_quality_60_unique.sam.out_FPKM.txt',row.names=T,sep='\t',quote=F)

#out <- merge(tpm_table, fpkm_table,
#                           by = 'row.names', all = TRUE)


outfile1=paste(c(basename(RC_file),"_TPM.txt"),collapse='')
outfile2=paste(c(basename(RC_file),"_FPKM.txt"),collapse='')
write.table(tpm_table, outfile1,row.names=F,sep='\t',quote=F)
write.table(fpkm_table, outfile2,row.names=F,sep='\t',quote=F)


#######################################
#analyze reps to get median/mean values
#######################################
# load packages
library(dplyr)
install.packages("matrixStats")
library(matrixStats)
# read in file with column name // tissue type
# NOTE: only use for samples that have rep- where reps have the same tissue type
tissue_file <- "tissue-rep_list.txt"
reps <- read.table(tissue_file,head=F,sep='\t',stringsAsFactors=F,row.names=NULL)
tissue_list<- unique(reps$V2) # get unique tissue values
# get median 
tbl_median_counts_tpm <- data.frame(row.names = tpm_table[,1]) # start newtable with gene row names
tbl_median_counts_fpkm <- data.frame(row.names = fpkm_table[,1])
# loop through tissues to get reps, then get median of reps and add to table
for (i in 1:length(tissue_list)) {
  subdat<- reps %>% filter(V2 == tissue_list[i]) #get list of reps for tissue type
  newdata1 <- tpm_table[,c(subdat$V1)] #subset reps
  newdata2 <- fpkm_table[,c(subdat$V1)]
  if(!is.null(ncol(subdat))){
  mat_num1 <- matrix(as.numeric(newdata1),    # Convert to numeric matrix
                  ncol = ncol(newdata1))
  mat_num2 <- matrix(as.numeric(newdata2),    # Convert to numeric matrix
                     ncol = ncol(newdata2))
  med_dat1 <- as.data.frame(rowMedians(mat_num1),optional=T) #get median - for mean can use rowMeans
  med_dat2 <- as.data.frame(rowMedians(mat_num2),optional=T)
  colnames(med_dat1)= paste(tissue_list[i],sep='') # add column name
  colnames(med_dat2)= paste(tissue_list[i],sep='')
  tbl_median_counts_tpm <- cbind(tbl_median_counts_tpm, med_dat1) # add median to table
  tbl_median_counts_fpkm <- cbind(tbl_median_counts_fpkm, med_dat2)
  } 
}

# run following line if you have sample without rep, replace peduncle with name and 4 with column index
tbl_median_counts_tpm <- cbind(peduncle=tpm_table[,4],tbl_median_counts_tpm)
tbl_median_counts_fpkm <- cbind(peduncle=fpkm_table[,4],tbl_median_counts_fpkm)

# write out median TPM
tbl_median_counts_tpm <- cbind(gene=tpm_table[,1],tbl_median_counts_tpm)
write.table(tbl_median_counts_tpm,'Median_TPM.txt',sep='\t',row.names=F,quote=F)
# write out median FPKM
tbl_median_counts_fpkm <- cbind(gene=fpkm_table[,1],tbl_median_counts_fpkm)
write.table(tbl_median_counts_fpkm,'Median_FPKM.txt',sep='\t',row.names=F,quote=F)
