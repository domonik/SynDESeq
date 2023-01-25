
suppressPackageStartupMessages({
  require(clusterProfiler)
  require(org.Ssp.eg.db)
})
print(keytypes(org.Ssp.eg.db))

defile <- snakemake@input[["defile"]]
ontology <- snakemake@wildcards[["ontology"]]

detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
#detable <- detable[!grepl("UTR", rownames(detable)), ]
#rownames(detable) <- gsub("_5UTR", "", rownames(detable))
#print(detable)
up <- detable$log2FoldChange >= 0.8 & detable$padj < 0.05
up <- detable[up, ]
down <- detable$log2FoldChange <= -0.8 & detable$padj < 0.05
down <- detable[down, ]
egoBPup <- enrichGO(gene = as.character(rownames(up)),
                    universe = as.character(rownames(detable)),
                    OrgDb = org.Ssp.eg.db,
                    keyType = "GID",
                    ont = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 100,
                    maxGSSize = 500,
                    readable = TRUE)
summary <- data.frame(egoBPup )
write.table(summary , file = snakemake@output[["up"]], row.names=FALSE, sep="\t")

egoBPdown <- enrichGO(gene = as.character(rownames(down)),
                    universe = as.character(rownames(detable)),
                    OrgDb = org.Ssp.eg.db,
                    keyType = "GID",
                    ont = ontology,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 100,
                    maxGSSize = 500,
                    readable = TRUE)

summary <- data.frame(egoBPdown )
write.table(summary , file = snakemake@output[["down"]], row.names=FALSE, sep="\t")