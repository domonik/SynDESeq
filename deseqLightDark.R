
suppressPackageStartupMessages({
  require(DESeq2)
  require(pheatmap)
})



countFile <- snakemake@input[["counts"]]
annotationFile <- snakemake@input[["annotation"]]
normGenesFile <- snakemake@input[["norm_genes"]]

normGenes <- read.table(normGenesFile, header=F, sep="\t", row.names=1, check.names=F)
countData <- as.matrix(read.table(countFile, header = T,  sep="\t", row.names=1, check.names = F, comment.char = ""))
annotationData <- as.matrix(read.table(annotationFile, header = T, row.names = 1, sep="\t", check.names = F, comment.char = ""))

formula <- as.formula(snakemake@params[["design"]])
dds <- DESeqDataSetFromMatrix(countData = countData, colData = annotationData, design = formula)

normGenes <- rownames(countData) %in% rownames(normGenes)

dds <- estimateSizeFactors(dds, controlGenes=normGenes)
dds <- DESeq(dds)

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