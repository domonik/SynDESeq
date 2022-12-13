suppressPackageStartupMessages({
    require(DESeq2)
    require(pheatmap)
    require(ggplot2)
})


countFile <- snakemake@input[["counts"]]
annotationFile <- snakemake@input[["annotation"]]
#countFile <- "Data/JoinedCounts.csv"
#annotationFile <- "Data/annotation.csv"

countData <- as.matrix(read.table(countFile, header = T,  sep="\t", row.names=1, check.names = F, comment.char = ""))
annotationData <- as.matrix(read.table(annotationFile, header = T, row.names = 1, sep="\t", check.names = F, comment.char = ""))
keep <- rowSums(countData) >= 10
countData <- countData[keep, ]
ercc <- grepl("ERCC", rownames(countData))
formula <- as.formula(snakemake@params[["design"]])
dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)

#dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = ~sample_group)
#dds <- DESeq(dds, controlGenes=ercc)
dds <- estimateSizeFactors(dds, controlGenes=ercc)
dds <- DESeq(dds)
print(sizeFactors(dds))
print(colSums(counts(dds)))
#dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind=TRUE)
pca_data <- plotPCA(rld, intgroup="sample_group", returnData=TRUE, ntop=500)
write.csv(pca_data, snakemake@output[["pca_data"]], sep="\t")
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
#write.csv(correlations, snakemake@output[["correlations"]])
pheatmap(rld_cor, filename=snakemake@output[["heatmap"]])
png(snakemake@output[["dispest"]])
plotDispEsts(dds)
dev.off()
save(dds, file=snakemake@output[["deseq_result"]])




