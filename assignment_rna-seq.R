# Elements of bioinformatics, RNA-seq data differential expression (DE) analysis
# In this assignment, we will compare the gene expression abundances in lung and lymph human samples.
# This script is adapted from DESeq2 method manual, and the utilized data is available on CSC Chipster (https://chipster.csc.fi/). 

###################### Problem 1 #######################
###########################################################

# Important: when unfamiliar with a function, use ? to search R documentation for help and check the help window.
# Example: ?head
# Use # sign to write your notes and your answers for problem 1.

setwd("/home/rna-seq/")

# Load/install required packages:
install.packages('tinytex')
tinytex::install_tinytex()
install.packages("pheatmap")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("apeglm") # This might take 1 to 2 minutes!

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(apeglm)


# Load table of counts (expression abundances for all experimental samples):
exp.levels <-  read.csv("counttab.csv", row.names = 1)

# Check how does the data looks like:
head(exp.levels)

# As you can see, rows represent genes, column 1 to column 4 represent gene position on the genome, 
# and columns 5 to 14 represent the samples and corresponding quantified expression abundances.

# How to access a specific gene expression level for all samples (how to access a row from a matrix/data frame in R)?
exp.levels["ENSG00000000971",]

# How to access the expression level for gene ENSG00000000971 from lymph_4 sample:
exp.levels["ENSG00000000971","lymph_4"]

# You can also access the chromosome for 10th gene as below:
exp.levels$chr[10]

# In order to run differential expression analysis, we need to remove columns with information about gene position, 
# meaning removing columns 1 to 4:
count.tab <-  exp.levels[,-c(1:4)] # From exp.levels data frame, keep all the rows but remove columns 1 to 4.

# Now you have the table of counts ready for DE analysis. 
# Check your dimension using dim() function. How many samples and genes are there?
dim(count.tab)
# Here, the count table includes 60675 genes and 10 samples (5 lung samples and 5 lymph samples).

# Load the annotation file:
annot.tab <-  read.csv("annottab.csv", row.names = 1)

# Check the annotation file:
annot.tab
# Its a metadata about samples available in our experiment. 
annot.tab$tissue = factor(annot.tab$tissue)

# Check if the annotation and count tables are consistent:
all(row.names(annot.tab) == colnames(count.tab))
# If you get TRUE, they are in the same order, otherwide you need to reorder the tables to match each other. 

# Create DESeq2 input object (DESeqDataSet):
deseqmat <-  DESeqDataSetFromMatrix(countData = count.tab,
                                  colData = annot.tab,
                                  design = ~ tissue)
# Here, the design formula is the  variables which will be used in modeling: tissue of interest in case of our analysis.

# For detailed information please check the package documentation:
?DESeqDataSetFromMatrix

# Check the DESeqDataSet:
deseqmat

# The DESeqDataSet is an R object to store the read counts and the intermediate results based on your analysis.
# counts() function can be used the expression levels from deseqmat object:
head(counts(deseqmat))

# Filter the genes with low expression abundances:
# As you saw, the count table includes 60675 genes out of which, many of them could have expression level equal or close to zero.
# In this step, we try to remove the genes with low expression level across all samples.
# This will help reducing the data matrix dimension and improving the analysis running time.
low.genes <-  rowSums(counts(deseqmat)) < 5
# Here, counts() function provides access to table of count from DESeqDataSet object.
# This way, we find the genes that include less than 5 reads for all experimental samples.

# Check how many of the genes have less than 5 reads in total:
table(low.genes)

deseqmat <-  deseqmat[!low.genes,]

# Check the dimension of count data matrix after filtering the lowly expressed genes:
dim(counts(deseqmat))

# Check the sample level quality and separation of samples based on experimental design:
# In order to do so, the data needs to be transformed (check the variance stablization using DESeq2) and normalized:
var.stablized <-  vst(deseqmat)

# Check the documentation for vst() function: ?vst

# Now, you can make the PCA plot using plotPCA() function available in DESeq2:
DESeq2::plotPCA(var.stablized, intgroup = "tissue")

# Based on the PCA plot, what do you think about the experimental design set? 

# Note: function plotPCA() is available in multiple packages with the same name (check ?plotPCA).
# To access this function from that specific package you aim (DESeq2 in this case), you need to provide the package name followed by 
# double colons before the function names. 

# Lets check how the raw and processed data differ from each other:
# raw counts:
head(counts(deseqmat))
# normalized and transformed counts:
head(assay(var.stablized))

# It is also possible to check the distribution of expression levels for each sample using boxplots:
# Samples distribution for raw data:
boxplot(log10(counts(deseqmat)), las = 2)
# Samples distribution for normalized and transformed data:
boxplot(assay(var.stablized), las = 2)
# These boxplots demonstrate the overall density distribution is identical after normalization and transformation for all samples. 

# Additionally, one can use hierarchical clustering techniques to check if the sample groups are separated from each other,
# or if any outlier exists:
# First calculate the distance matrix between samples using dist() function:
sample.distances <-  dist(t(assay(var.stablized)))
# you can check the dist() function documentation: ?dist

# Check how the distance matrix looks like:
sample.distances

# Now that we have the distances matrix we can make a heatmap visualization to present it.
# Here, we will use pheatmap() from pheatmap package for this purpose:
colors <-  colorRampPalette(brewer.pal(9, "PuBu"))(100)
pheatmap(as.matrix(sample.distances),
         clustering_distance_rows = sample.distances,
         clustering_distance_cols = sample.distances,
         col = colors)
# You can check the heatmap() function and arguments from the package documentation: ?pheatmap

# What do you think about the experimental design of this study? 
# What is the main source of variance among samples?

# Now that we explored the data and quality checks look promising, it is time to perform the actual testing to 
# detect the differentially expressed genes: 

# Run differential expression testing using DESeq2:
# In DESeq2, a single function (DESeq()) is responsible for performing normalization, variance stablization, model fitting and testing.
# Of course the are also separate functions to do these tasks as well.
de.test <-  DESeq(deseqmat)
de.test

# Obtain the DESeq2 result matrix using the results() function:
res <-  results(de.test)
# Check the results() function and its arguments for more details: ?results
# Lets check how the result matrix looks like:
head(res)

# Lets check the genes with pvalue < 0.01:
res.pval01 <-  subset(res, pvalue < 0.01)
# subset() function returns subsets of matrix or data frames that fulfill the determined criteria (here pvalue < 0.01)
summary(res.pval01)
# Can you check how many of the genes are determined are significant with pvalue < 0.05? 
# How many of them are up-regulated?

# Can you check how many of the genes have adjusted pvalue less than 0.05 using subset() function?
# Your answer here:


# Visualizing the expression levels for interesting (significant) genes:
# Lets suppose we are interested to investigate the expression level of a single gene across the study groups.
# plotCounts() function shows the normalized and transformed (log) expression levels of a specific gene:
plotCounts(de.test, gene=which.min(res$padj), intgroup="tissue")
# Here we plotted the normalized counts for gene with smallest fdr on log scale.
plotCounts(de.test, gene=which.max(res$log2FoldChange), intgroup="tissue")
# please check the arguments for plotCounts() function for more details.

# The plotCounts() function can shows the expression levels only for one gene.
# If you are interested to check the expression levels of a set of genes across sample groups, heatmaps are one suitable option.

# Lets assume you are interested to check the expression levels for 20 first genes with highest positive log fold change.
# We can check this using heatmaps.
# First lets find the target genes by sortung the DESeq2 results matrix based on the log fold change:
ordered.res <-  res[order(res$log2FoldChange, decreasing = TRUE),]
# Then select the first 20 genes:
top.genes <-  row.names(ordered.res)[1:20]

# Now we can plot the heatmap for the selected genes using the normalized and transformed expression levels:
pheatmap(assay(var.stablized)[top.genes,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE)

# Based on the heatmap, in which sample group the selected genes are over-expressed?

# More exploration of the results: MA plot
# plotMA() function outputs a scatter plot of log fold changes versus the mean of normalized expression values.
# The gray dots are the non-significant genes and the blue dots are the significantly differentially expressed genes (default fdr = 0.1).
plotMA(res, ylim=c(-3,3))

# We discussed that DESeq2 aims to shrink the log fold changes for very lowly expressed genes.
# Here, using lfcShrink() function, we can calculate the shrunken log fold changes and plot them versus the mean expression levels again:
# Check the MA plot for the shrinked log fold changes:
resLFC <-  lfcShrink(de.test, coef="tissue_lymph_vs_lung")
plotMA(resLFC, ylim=c(-3,3))
# Do you know why DESeq2 aims to shrink the log fold changes?

# Now you can save your analysis results (all R objects in current session) using save.image() command:
save.image(file = "/home/rstudio/rna-seq-analysis.RData")
# You can download this image from the file browser tab (the lower right window) by selecting the file -> "More" -> "Export"
# Note the file is located at /home/rstudio/.

# You can write the DESeq2 results matrix to a file for further use as follow:
write.csv(res, file = "/home/rstudio/deseq2_de_results.csv")
# You can downlad this file the same as the R image.

# To save your script: File (from the menu bar) -> Save as -> File name: /home/rstudio/assignment
# You can download this script from the file browser tab (the lower right window) by selecting the file -> "More" -> "Export"

###########################################################
###########################################################

###################### Problem 2 #######################
###########################################################

# A) Can you check the sample distances using "manhattan" method and plot the hierarchical clustering of samples again?
# Your code here:


# B) Can you investigate the MA-plots to check for log fold change versus mean expression levels and this time
# set the adjusted pvalue threshold to 0.05?
# You can check the plots for both raw and shrunken fold changes.
# Your code here:



# C) Please define a desirable criteria to check the expression level of a single gene across sample groups.
# For example, you can find a gene with largest minus fold change and show its expression level across tissues using plotCounts() function.
# Your code here:



# D) Please define a desirable criteria to check the expression level of a set of genes across sample groups.
# For example, you can select top 20 genes with smallest adjusted pvalues and use pheatmap function to check their expression levels across sample groups.
# Your code here:




# Use same commands save(), save.image() to save your analysis results and same approach as in problem 1 to downlowd them to your computer.


# To save your script: File (from the menu bar) -> Save as -> File name: /home/rstudio/assignment
# You can download this script from the file browser tab (the lower right window) by selecting the file -> "More" -> "Export"

###########################################################
###########################################################

# Note when the analysis is finished you can compile a report using: File -> Compile Report... -> choose format of interest -> Compile 
# Please remember to save the report on your computer.
