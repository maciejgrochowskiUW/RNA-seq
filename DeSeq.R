library ("DESeq2")
setwd("~/Documents/workdir/mapped/htseq")
#first u need to point the directory with your htseq-count output data
directory <- "/home/a/Documents/serv/htseq"


#specify which files to read in using list.files, and select those files which contain the string "counttable" using grep.
#The sub function is used to chop up the sample filename to obtain the condition status, alternatively you might read in a phenotypic table using read.table.

sampleFiles <- grep("counttable",list.files(directory),value=TRUE)
sample <- sub("counttable_","",sampleFiles)
sample2 <- sub(".txt","",sample)
sampleCondition <- sub("_[[:digit:]].txt","",sample)
sampleTable <- data.frame(type = sample2,
                          fileName = sampleFiles,
                          condition = sampleCondition)
rownames(sampleTable) <- sampleTable$type
#now u can load your data to R with DESeqDataSetFromHTSeqCount command
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
#add aditional column with sample name
colData(ddsHTSeq)[,2] <- rownames(colData(ddsHTSeq))
names(colData(ddsHTSeq))[2] <- "type"
colData(ddsHTSeq)$type <- as.factor(colData(ddsHTSeq)$type)


ddsHTSeq
#it should generate something like this
#class: DESeqDataSet 
#dim: 5137 4 
#metadata(1): version
#assays(1): counts
#rownames(5137): SPAC1002.01 SPAC1002.02 ... SPMTR.03 SPMTR.04
#rowData names(0):
#colnames(4): dis3l2_S16 dis3l2_S17 wt_S14 wt_S15
#colData names(2): condition type

#Get rid of genes that have less than 20 reads in all samples
keep <- rowSums(counts(ddsHTSeq)) >= 20
dds <- ddsHTSeq[keep,]
dds

#set your reference sample
dds$condition <- relevel(dds$condition, ref = "WT")

#DESeq2 provides a function collapseReplicates which can assist in combining the counts
#from technical replicates into single columns of the count matrix.
#The term technical replicate implies multiple sequencing runs of the same library.
#You should not collapse biological replicates using this function.

#The standard differential expression analysis steps are wrapped into a single function, DESeq.
dds <- DESeq(dds)
res <- results(dds)
res

#Log fold change shrinkage for visualization and ranking. Use command below
resultsNames(dds)
#and it will print something like this
#[1] "Intercept"              "condition_dis3l2_vs_wt". Now use this output for next command
resLFC <- lfcShrink(dds, coef="condition_C1_vs_WT", type="apeglm")
resLFC

#order results table by the smallest p value
resOrdered <- res[order(resLFC$pvalue),]
summary(res)

#How many adjusted p-values were less than 0.1
sum(res$padj < 0.1, na.rm=TRUE)

#The results function contains a number of arguments to customize the results table which is generated.
#If the adjusted p value cutoff will be a value other than 0.1, alpha should be set to that value.
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("type"))
plotPCA(vsd, intgroup=c("condition"))


plotMA(res, ylim=c(-3,3))

#It is more useful to visualize the MA-plot for the shrunken log2 fold changes,
#which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.
plotMA(resLFC, ylim=c(-3,3))

##########DOESNT WORK#####################
#library("pheatmap")
#ntd <- normTransform(dds)
#select <- order(rowMeans(counts(dds,normalized=TRUE)),
#               decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("condition","sample")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)

#library("IHW")
#resIHW <- results(dds, filterFun=ihw)
#summary(resIHW)
#sum(resIHW$padj < 0.1, na.rm=TRUE)
#metadata(resIHW)$ihwResult




#write.csv(as.data.frame(resOrdered), 
#         file="deseq2_dis3l2_results.csv")

df1 <- as.data.frame(resOrdered)
DIS3L2_mysubset <- subset(df1, (log2FoldChange > 1 | log2FoldChange < -1) & padj < 0.05)
DIS3L2_mysubset2 <- subset(DIS3L2_mysubset, (log2FoldChange > 1))
write.csv(DIS3L2_mysubset2, file = "deseq_DIS3L2_up.csv")
DIS3L2_mysubset3 <- subset(DIS3L2_mysubset, (log2FoldChange < -1))
write.csv(DIS3L2_mysubset3, file = "deseq_DIS3L2_down.csv")




d <- plotCounts(dds, gene="SPAC2C4.07c", intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0))

coldata <- data.frame(x=factor(c(1,1,2,2)))
cts <- matrix(1:16,ncol=4)
dds2 <- DESeqDataSetFromMatrix(cts, coldata, ~x)
sizeFactors(dds2) <- rep(1,4)
plotCounts(dds2, 1, "x")
abline(h=counts(dds2, normalized=TRUE)[1,] + .5)