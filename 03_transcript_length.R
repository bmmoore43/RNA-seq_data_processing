if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

library("GenomicFeatures")

txdb <- makeTxDbFromGFF("Jascendensv1.1.gene.gff3",format="gff")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
cat("gene\ttranscript_length\n", file="Jascendensv1.1.gene.gff3.transcript_length.txt")
write.table(exonic.gene.sizes[2,], file="Jascendensv1.1.gene.gff3.transcript_length.txt",row.names = TRUE, sep="\t", append=TRUE,quote=FALSE,)
