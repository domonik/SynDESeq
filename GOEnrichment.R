
suppressPackageStartupMessages({
  require(clusterProfiler)
  require(org.Ssp.eg.db)
})
print(keytypes(org.Ssp.eg.db))

defile <- snakemake@input[["defile"]]

detable <- read.table(defile,sep="\t",header=TRUE)
detable <- na.omit(detable)
#detable <- detable[!grepl("UTR", rownames(detable)), ]
#rownames(detable) <- gsub("_5UTR", "", rownames(detable))
#print(detable)
up <- detable$log2FoldChange >= 1 & detable$padj < 0.05
up <- detable[up, ]
down <- detable$log2FoldChange <= -1 & detable$padj < 0.05
down <- detable[down, ]
egoBPup <- enrichGO(gene = as.character(rownames(up)),
                    universe = as.character(rownames(detable)),
                    OrgDb = org.Ssp.eg.db,
                    keyType = "GID",
                    ont = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 100,
                    maxGSSize = 500,
                    readable = TRUE)
print(egoBPup)
print(egoBPup$Description)
