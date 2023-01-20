if (!require("BiocManager", quietly = TRUE)){

    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

}
if (!require('tidyr', quietly = TRUE)){

    install.packages('tidyr')
}

BiocManager::install(c("DESeq2", "pheatmap"))
BiocManager::install(c("clusterProfiler" ))
BiocManager::install(c("AnnotationForge" ))

file.create(snakemake@output[["install_file"]])
