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

#RC_file <- "HTSeqCount_9-month-old-internode-2_sam3_sorted_quality_60_unique.sam.out"
#TL_file <- "Jascendensv1.1.gene.gff3.transcript_length.txt"

Exp <- read.table(RC_file,head=F,sep='\t',stringsAsFactors=F,row.names=1)
dat <- read.table(TL_file,head=T,sep='\t',stringsAsFactors=F,row.names=1)
#Expression <- cbind(Exp[order(row.names(Exp)),],'Length'=dat[order(row.names(dat)),])
#Expression <- as.data.frame(Expression)
Expression <- merge(Exp, dat,
                          by = 'row.names', all = TRUE)
Expression <- Expression[6:29125,]
#samp2 <- Expression[,-1]
#rownames(samp2) <- Expression[,1]
Len <- Expression$transcript_length
tpm_table <- calculateTPM(as.matrix(Expression[,2:(ncol(Expression)-1)]), Len)
rownames(tpm_table) <- Expression[,1]
#write.table(tpm_table,'HTSeqCount_9-month-old-internode-2_sam3_sorted_quality_60_unique.sam.out_TPM.txt',row.names=T,sep='\t',quote=F)

fpkm_table <- calculateFPKM(as.matrix(Expression[,2:(ncol(Expression)-1)]), Len)
rownames(fpkm_table) <- Expression[,1]
#write.table(fpkm_table,'HTSeqCount_9-month-old-internode-2_sam3_sorted_quality_60_unique.sam.out_FPKM.txt',row.names=T,sep='\t',quote=F)

out <- merge(tpm_table, fpkm_table,
                           by = 'row.names', all = TRUE)

colnames(out)<- c("gene","TPM","FPKM")
outfile=paste(c(basename(RC_file),"_TPM-FPKM.txt"),collapse='')
write.table(out, outfile,row.names=F,sep='\t',quote=F)


###
# tissue_list <- unique(str_remove_all(as.character(colnames(Expression)), "_\\d$"))
# tbl_median_counts <- data.frame(row.names = rownames(tpm_table))
# for (i in 1:length(tissue_list)) {
#   subdat <- tpm_table[,grep(colnames(tpm_table),pattern=tissue_list[i],fixed = TRUE)]
#   if(!is.null(nrow(subdat))){
# 	  new_dat <- as.data.frame(rowMedians(subdat), optional = T)
# 	  colnames(new_dat) = paste(tissue_list[i],sep='')
# 	  tbl_median_counts <- cbind(tbl_median_counts, new_dat)
# 	}
# }
# write.table(tbl_median_counts,'Median_TPM.txt',sep='\t',quote=F)
