library(AnnotationForge)
library("tidyr")


mapping_data <- read.table(snakemake@input[["locustag2GO"]], sep="\t", header=TRUE)
mapping_data$EVIDENCE <- "experimental"
colnames(mapping_data) <- c("GID", "GO", "EVIDENCE")
symbol_data <- read.table(snakemake@input[["locustag2Symbol"]], sep=",", header=TRUE, fill=TRUE)
colnames(symbol_data) <- c("GID", "SYMBOL")
symbol_data <- symbol_data %>% drop_na()
dir.create(snakemake@output[["annodb"]], recursive=TRUE)
rv <- makeOrgPackage(gene_info=symbol_data, go=mapping_data, version="0.1", maintainer="DR <schdruzzi@gmail.com>",
                     author="DR <schdruzzi@gmail.com>",
                     outputDir = snakemake@output[["annodb"]],
                     tax_id = "42", genus = "Synechocystis", species="sp", goTable = "go")
install.packages(rv, repos=NULL)
file.create(snakemake@output[["finished_file"]])





