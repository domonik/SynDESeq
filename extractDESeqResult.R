

suppressPackageStartupMessages({
  require(DESeq2)
})


condition <- snakemake@wildcards[["condition"]]
baseline <- snakemake@wildcards[["baseline"]]
file <- snakemake@input[["deseq_result"]]
variable <- snakemake@params[["design"]]
print(file)
load(file=file)
contrast <- c(variable, condition, baseline)
result_table <- results(dds, contrast=contrast, alpha = 0.05)
write.table(result_table, snakemake@output[["result_table"]], sep="\t")