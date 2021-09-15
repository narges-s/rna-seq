# Loading libraries:
library(DESeq2)
library(airway)
library(pasilla)


# DESeq2 package:
# airway package:
# pasilla package:
# We analyze RNA-seq data from x paper. The data is available as a bioconductor package.

# DESeq2 expects data from ... to be like ...


# Load the input data and make some data preparation:

count.table <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
annotation <- system.file("extdata",
                       "pasilla_sample_annotation.csv",
                       package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(count.table,sep="\t",row.names="gene_id"))
coldata <- read.csv(annotation, row.names=1)
coldata <- coldata[,c("condition","type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)



# Check the count table (gene expression abundances):
head(cts)


# Check the meta data matrix:
coldata

# Explain about the meta data matrix:


# Make the namaing and colnames of count table consistent:
row.names(coldata) = sub("fb", "", rownames(coldata))


all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

dds = DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds


dds <- DESeq(dds)
res <- results(dds)
res

# res <- results(dds, name="condition_untreated_vs_treated")
res <- results(dds, contrast=c("condition","treated","untreated"))

resultsNames(dds)
# [1] "Intercept"                      "condition_untreated_vs_treated"


sum(res$pvalue < 0.01, na.rm=TRUE)
# [1] 581

pdf("maplot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()


res05 <- results(dds, alpha=0.05)
summary(res05)


resNorm <- lfcShrink(dds, coef=2, type="normal")

pdf("test_1.pdf")
plotMA(resNorm, main="normal")
dev.off()

pdf("test_1.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
dev.off()

write.csv(as.data.frame(res05), 
          file="condition_treated_results.csv")

library("pheatmap")

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])

pdf("test2.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])

pdf("test2.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()


ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])

pdf("test2.pdf")
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()


vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("test_3.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf("test.pdf")
plotPCA(vsd, intgroup=c("condition", "type"))
dev.off()

library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf("test.pdf")
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
dev.off()


par(mar=c(8,5,2,2))
pdf("test.pdf")
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()

pdf("test.pdf")
plotDispEsts(dds)
dev.off()


use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
#Histogram of p values for all tests. The area shaded in blue indicates the subset of those that pass the filtering, the area in khaki those that do not pass:

pdf("test.pdf")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))
dev.off()


# *********************************************************************************** #
# Install if not available:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("pasilla")
# BiocManager::install("airway")
# library(airway)
# data(airway)
# airway$SampleName
# airway$cell
# airway$Run
# airway$Sample
# temp = airway@rowRanges
# data(gse)


# Lecture slides figures:
# Sample distances plot:
require(airway)
data(gse)
gse$cell <- gse$donor
gse$dex <- gse$condition
levels(gse$dex) <- c("untrt", "trt")
gse$dex <- relevel(gse$dex, "untrt")

library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ cell + dex)

countdata <- round(assays(gse)[["counts"]])

coldata <- colData(gse)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ cell + dex)

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)
sampleDists <- dist(t(assay(vsd)))

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
samplenames = vsd$dex
# [1] untrt trt   untrt trt   untrt trt   untrt trt
samplenames = c("WT_1","MUT_1","WT_2","MUT_2","WT_3", "MUT_3","WT_4", "MUT_4")

rownames(sampleDistMatrix) <- samplenames
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pdf("sample.dist.pdf")
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()



# MDS plot:
# If really necessary




