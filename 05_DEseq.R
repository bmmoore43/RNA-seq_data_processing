###################################
# Differential expression (edgeR) #
###################################
# install edgeR
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

# load the package and data
# Run every time to you start a new R session.
library(edgeR)

# Download the expression data  
counts <- read.table('counts_data.csv', header=T, row.names='gene', sep='\t')
head(counts)
#summary(counts$control_1) # Get a summary of the expression levels for the first replicate of the control treatment

# You first need to tell edgeR where to look for your data (i.e. the genes, counts, and treatment groups) using DGEList
# DGEList is a data container used by edgeR to organize this information
# To define the treatment groups, you need to tell EdgeR what columns belong to the same treatment group (meaning they are replicates) 

treatments <- c("C","C","C","Cold","Cold","Cold") # you are making a new list called "group" using c() of your treatments
# I think this can be read in as a table
d <- DGEList(counts, group=treatments, genes=rownames(counts)) # rownames(counts) tells DGEList that the gene names are rownames in your counts dataframe

# d is now the variable that contains your DGEList
# You can visualize the frequency distribution of your counts in your samples using a histogram.
# Due to the large variability in counts per gene it is useful to use the log transform counts for normalization.
hist(d$counts, breaks = 50, xlab = "Counts")
hist(log2(d$counts), xlim=c(0, 20), breaks=50, xlab = "LogCounts")


# Normalize between samples
# The goal is to normalize the counts across samples so that MOST genes have the same abundance across the different samples. 
# The idea is that while you expect some genes will be differentially expressed across the datasets, most genes will not be affected
# by your treatment. We can use those "unchanging" genes to normalize between the samples. 
d <- calcNormFactors(d) #This calculates the effective library size and normalization factors. 
d$samples 



# Another way to look at the variation between your samples is to look at the "Biological Coefficient of Variation" (BCV)
# Because of the way variance is calculated, genes with higher expression would have higher variance.
# The BCV corrects for this by taking mean expression level into consideration. 
# The BCV is determined by estimating the common (across all genes in all samples) and 
# tagwise (across all samples for an individual gene) dispersal. 
# A common BCV between 0.2-0.4 is considered good enough to be able to identify differentially expressed genes. 

# Estimate the dispersion (i.e. the BCV between biological replicates).
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)


# Counts per million (CPM) mapped reads are counts scaled by the number of fragments you sequenced (N) times one million.
# CPM is not scaled by the gene length. 
# Average log CPM : average log CPM in each sample for each gene.



###### Use the Fisher's Exact test to look for differentially expressed genes. #######

# First we'll look for genes differentially expressed between control & cold treatment

# The Fisher's exact test takes into account the "Tagwise Dispersian" when calculating the p-value
c_vs_cold <- exactTest(d, pair=c('C','Cold')) 

# Since we are calculating significance for thousands of genes, we need to use false discovery rate correction
# Here we use the Benjamini-Hochberg (BH) method (BH) to correct the p-values, the correct p-values are called q-values.
# We are saving the corrected p-values in a new column in the c_vs_cold$table called "FDR"
c_vs_cold$table$FDR <- p.adjust(c_vs_cold$table$PValue, method='BH') 

# Select the genes that are significantly differentially expressed (q-value <= 0.05) and have a logFC > 2 or < -2
c_vs_cold_genes <- subset(c_vs_cold$table, c_vs_cold$table$FDR <= 0.05 & abs(c_vs_cold$table$logFC) >=2)
write.table(data.frame("gene"=rownames(c_vs_cold_genes),c_vs_cold_genes),file="Control_vs_Cold_DEgenes.txt",quote=FALSE,sep="\t",row.names=FALSE)


# Check out the relationship between logFC and significance by making a "volcano" plot:
plot(-log(c_vs_cold_genes$FDR)~c_vs_cold_genes$logFC, 
     ylab= "logFDR",
     xlab = "logFC") 

# Save your results (you do not need to turn these in, but you may want to refer to these results later)
write.table(data.frame("gene"=rownames(c_vs_cold$table),c_vs_cold$table),file="Control_vs_Cold.txt",quote=FALSE,sep="\t",row.names=FALSE) #quote=FALSE, to remove the "" that are automatically added on. 

