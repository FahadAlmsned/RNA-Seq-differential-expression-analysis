library(DESeq2)
library(IHW)
library(dplyr)
library(devtools)
library(rafalib)
library(vsn)
library(hexbin)
library(Rsamtools)
library(annotables)
library(biobroom)
library(ggplot2)
library(pheatmap)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(pathview)
library(gage)
library(gageData)
#--------------

counts = read.delim("counts.txt",sep = "")
sample.table = read.csv("samples.csv",header=TRUE, row.names= 1)

#Experimental design and running DESeq2
dds <- DESeqDataSetFromMatrix(counts, colData=sample.table, design= ~Dis)

#filtering
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]

# Size factors and colsum
sizeFactors(dds)
colSums(counts(dds))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#Levels
design(dds)
dds$Dis <- relevel(dds$Dis, ref = "N")
levels(dds$Dis)

#runs the DESeq2 model
dds <- DESeq(dds)
res <- results(dds)
head(res)
summary(res)
sum(res$padj < 0.1, na.rm = TRUE)
table(res$padj < 0.1) #Adjusted p-value


#More information on results columns
mcols(res)$description

#ordr the results
esorder = res[order(res$padj),]
head(resorder)

#Exporting results to CSV files
write.csv(as.data.frame(resorder),file = "MS_N.csv")


# alpha Independent filtering is further 
res05 <- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE)


#plotMA shows the log2 fold changes  over the mean of normalized counts for all the samples in the DESeqDataSet.
plotMA(res, ylim = c(-10,10))

#interactively detect the row number of individual genes by clicking on the plot
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

#moderation of log2 fold changes
resultsNames(dds)
resLFC <- lfcShrink(dds, contrast = c("Dis","MS","N"), res = res)
resLFC

# visualize the MA-plot for the shrunken log2 fold changes
plotMA(resLFC, ylim = c(-2,2))

#We can also test against a different null hypothesis.
res.thr <- results(dds, lfcThreshold = 0.2)
plotMA(res.thr, ylim = c(-10,10))

# volcano plot red if padj<0.05, orange of log2FC>1, green if both
with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "Volcano plot", xlim = c(-10,10)))
with(subset(res, padj < .1 ), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
with(subset(res, abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "orange"))
with(subset(res, padj < .05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "green"))

#A p-value histogram
hist(res$pvalue[res$baseMean > 1], col = "grey", border = "white", xlab = "", ylab = "", main = "")


#Examine the counts for the top gene, sorting by p-value:
plotCounts(dds, gene = which.min(res$padj), intgroup = "Dis")

#customized plotting
d <- plotCounts(dds, gene = which.min(res$padj), intgroup = "Dis",  returnData = TRUE)

ggplot(d, aes(x = Dis, y = count)) + geom_point(position = position_jitter(w = 0.1,h = 0)) + scale_y_log10(breaks = c(25,100,400))


#annotables

#Look at the human genes and trancrpits tables
grch38
grch38_tx2gene

grch38 %>% 
  dplyr::filter(biotype == "protein_coding" & chr == "1") %>% 
  dplyr::select(ensgene, symbol, chr, start, end, description) %>% 
  head

# tidy results with biobroom
res_tidy <- tidy.DESeqResults(res)
head(res_tidy)

n = res_tidy %>% 
  dplyr::arrange(p.adjusted) %>% 
  head(100) %>% 
  dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
  dplyr::select(gene, estimate, p.adjusted, symbol, description)

View(n)

#or out as text
#write.table(n, file="top_genes_MS_vs_N.txt")
write.csv(as.data.frame(n), file = "top_genes_MS_vs_N.csv")


#GAGE

columns(org.Hs.eg.db)

res$symbol = mapIds(org.Hs.eg.db,
                    keys = row.names(res), 
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")
res$entrez = mapIds(org.Hs.eg.db,
                    keys = row.names(res), 
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")
res$name =   mapIds(org.Hs.eg.db,
                    keys = row.names(res), 
                    column = "GENENAME",
                    keytype = "ENSEMBL",
                    multiVals = "first")

head(res, 10)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs, 3)
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# Get the results
keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)

# Look at both up (greater), down (less), and statatistics.
keg = lapply(keggres, head)
keg



# Get the pathways (change for less)
keggrespathways = data.frame(id = rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number() <= 5) %>% 
  .$id %>% 
  as.character()
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start = 1, stop = 8)
keggresids

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data = foldchanges, pathway.id = pid, species = "hsa", new.signature = FALSE)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data = foldchanges, pathway.id = pid, species = "hsa"))


#gene ontology
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchanges, gsets = gobpsets, same.dir = TRUE)

lapply(gobpres, head)
#


#In order to test for differential expression, we operate on raw counts and use discrete distributions as described in the previous section on differential expression.
#However for other downstream analyses - e.g. for visualization or clustering - it might be useful to work with transformed versions of the count data.


# Examine the log counts and the log normalized counts (plus a pseudocount)
log.norm.counts <- log2(counts(dds, normalized = TRUE) + 1) #log normalized counts (plus a pseudocount) 
#or
log.norm <- normTransform(dds) #metadata

rs <- rowSums(counts(dds))
mypar(1,2)
boxplot(log2(counts(dds)[rs > 0,] + 1))  # not normalized
boxplot(log.norm.counts[rs > 0,]) # normalized


#Extracting transformed values (better methods)
rld <- rlog(dds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
#head(assay(rld), 3)

#plotting
plot(log.norm.counts[,3:4], cex = .1)
plot(assay(rld)[,3:4], cex = .1)
plot(assay(vsd)[,3:4], cex = .1)

#Effects of transformations on the variance (standard deviation of rows over the mean)
meanSdPlot(assay(log.norm))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd))





#Data quality assessment by sample clustering and visualization

#Heatmap of the count matrix

#select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("Dis","Age","Ethnic")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)


#hierarchical clustering based on Euclidean distance matrix
par(mfrow=c(2,2))
plot(hclust(dist(t(assay(dds)))), labels=colData(dds)$Dis)
plot(hclust(dist(t(log.norm.counts))), labels=colData(dds)$Dis)
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$Dis)
plot(hclust(dist(t(assay(vsd)))), labels=colData(vsd)$Dis)


#(PCA) plot is a useful diagnostic for examining relationships
plotPCA(DESeqTransform(dds),intgroup=c("Dis"))
plotPCA(log.norm, intgroup=c("Dis"))
plotPCA(rld, intgroup=c("Dis"))
plotPCA(vsd, intgroup=c("Dis"))




