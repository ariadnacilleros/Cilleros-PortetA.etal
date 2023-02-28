### Classify imQTLs according Kim Hellmuth et al.
### Ariadna Cilleros-Portet

setwd("./")
colnames_vector <- c("phenotype_id","variant_id","tss_distance","af","ma_samples","ma_count","pval_g","b_g","b_g_se","pval_i","b_i","b_i_se","pval_gi","b_gi","b_gi_se")
library(data.table)

# Load imQTLs  
all <- fread("nominal_db_filterpval.txt", col.names = colnames_vector)
all$mqtl <- paste(all$phenotype_id, all$variant_id, sep="-")

# Load high proportion
high <- fread("highpercentile/nominal_db_grep.txt", col.names = colnames_vector)
high$mqtl <- paste(high$phenotype_id, high$variant_id, sep="-")

# Load low proportion
low <- fread("lowpercentile/nominal_db_grep.txt", col.names =  colnames_vector)
low$mqtl <- paste(low$phenotype_id, low$variant_id, sep="-")

# Get common imQTLs
common <- intersect(intersect(high$mqtl, low$mqtl), all$mqtl)
# 1,275 mQTLs in common

# Merge dfs
merged <- merge(x = high[,c("mqtl","b_g")], y = low[,c("mqtl","b_g")], by ="mqtl")
colnames(merged) <- c("mqtl","b_g_high","b_g_low")
merged_all <- merge(x = all[,c("mqtl","pval_gi","b_g", "b_gi")], y = merged, by ="mqtl")
rm(merged)

# Count positive, negatives and uncertain imQTLs 
uncertain <- merged_all[b_g_high > 0 & b_g_low < 0 | b_g_high < 0 & b_g_low > 0, ]
# 230 uncertain mQTLs

# Exclude uncertain imQTLs
merged_all <- merged_all[!(mqtl %in% uncertain$mqtl), ]

# Get positive imQTLs
positives <- merged_all[b_g_high > b_g_low,]
# 698 positives

# Get negative imQTLs
negatives <- merged_all[b_g_high < b_g_low,]
# 347 negatives


