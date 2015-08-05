#!/usr/bin/R

library(edgeR)

##Read in a single data file..  Each row is a transcript.  First column is IDs, each other column is raw counts from 1 sample.

rawdata <- read.delim("count.all.6.tissues")

dge <- DGEList(counts=rawdata[,c(2:52)], genes=rawdata[,1])

#Specify groups that each sample belongs to.

dge$samples$group <- c(rep("small intestine", 8), rep("heart", 9), rep( "prostate", 7), rep("liver", 5), rep("lymphnode", 13), rep("thyroid", 9))

dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

#Perform exact test, specifying pairwise comparison to make.

et <- exactTest(dge, pair=c("heart", "thyroid"))

#Adjust p value (BH method) and filter DE genes to <0.1.

signifAdjPValue = sum(p.adjust(et$table$table$PValue,method="BH") < 0.1)

signifList = topTags(et, n= signifAdjPValue)

#Sink output

sink("topPairwiseDEGenes.csv")
signifList
sink()


