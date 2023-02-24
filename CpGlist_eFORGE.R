### Get CpG list for eFORGE analysis
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# TensorQTL results
nominal <- fread("nominalfastqtl_format_pval5e-8.txt", col.names = c("cpg","sno","dist","p","b"))

# Filter CpGs with best p-value
library(dplyr)
nominal_cpg <- group_by(nominal, cpg)
nominal_cpg_filter <- filter(nominal_cpg, rank(p, ties.method="min")==1)

# Order by p-value column 
nominal_cpg_filter <- nominal_cpg_filter[order(nominal_cpg_filter$p), ]

# Write top 10,000
write.table(x = nominal_cpg_filter[1:10000,c("cpg","p")], file = "./cpg_list_eFORGE_top10000_sortedTOPpvalue.txt", quote = F, sep = "\t", row.names = F, col.names = F)
