### Get FastQTL format for mQTL database
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# Get TensorQTL results filtered by p-value < 5e-8
tensor <- fread("./nominal_db_filterpval.txt", 
                col.names = c("phenotype_id","variant_id","tss_distance","af","ma_samples","ma_count","pval_nominal","slope","slope_se"))


# Select columns according FastQTL format
tensor <- tensor[,c("phenotype_id","variant_id","tss_distance","pval_nominal","slope")]

# Write output file
fwrite(x = tensor, sep = "\t", quote = F, row.names = F, col.names = F, 
       file = "./nominal_db_filterpval_fastqtl.txt")
