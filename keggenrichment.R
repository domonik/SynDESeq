suppressPackageStartupMessages({
  require(clusterProfiler)
})


defile <- snakemake@input[["defile"]]

detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
detable <- detable[!grepl("UTR", rownames(detable)), ]

up <- detable$log2FoldChange >= 0.8 & detable$padj < 0.05
up <- detable[up, ]
down <- detable$log2FoldChange <= -0.8 & detable$padj < 0.05
down <- detable[down, ]
ekeggup <- enrichKEGG(gene         = as.character(rownames(up)),
                        universe      = as.character(rownames(detable)),
                        organism     = 'syn',
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
summary <- data.frame(ekeggup )
write.table(summary , file = snakemake@output[["up"]], row.names=FALSE, sep="\t")

ekeggdown <- enrichKEGG(gene         = as.character(rownames(down)),
                        universe      = as.character(rownames(detable)),
                        organism     = 'syn',
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)
summary <- data.frame(ekeggdown )
write.table(summary , file = snakemake@output[["down"]], row.names=FALSE, sep="\t")