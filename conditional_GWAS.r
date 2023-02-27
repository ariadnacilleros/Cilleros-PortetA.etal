### Conditional signals on GWAS
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# Load coloc results
coloc <- as.data.table(readRDS(file = "./output_coloc.RDS"))

### Subset coloc results by PPA4 > 0.8 (ref. Leyden et al.2022)
coloc <- coloc[coloc$PP.H4.abf > 0.8, ]
#2,177 coloc hits; 2,057 CpGs, 188 regions
fwrite(x = coloc, file = "./resultscoloc_SCZ.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# Load SMR+HEIDI results
smr <- fread("./nominal_inma_outcomeSCZ_multi.msmr")

### Subset SMR+HEIDI results 
smr$bonf <- p.adjust(smr$p_SMR_multi, "bonferroni")
smr_heidi_bonf <- smr[bonf < 0.05, ]
smr_heidi_bonf <- smr_heidi_bonf[p_HEIDI > 0.05 | is.na(p_HEIDI), ]

# Intersect SMR+HEIDI results and coloc 
intersect(smr_heidi_bonf$probeID, coloc$trait2)
# 96 CpGs/2,057 CpGs from colocalization

# Get COJO hits
smr_noheidi_bonf <- smr[bonf < 0.05, ]
smr_noheidi_bonf <- smr_noheidi_bonf[p_HEIDI < 0.05, ]
# 469 SMR hits with heterogeneity (p-HEIDI < 0.05)

# Intersect SMR hits with heterogeneity with coloc to perform COJO
cpgs_cojo <- intersect(smr_noheidi_bonf$probeID, coloc$trait2)

### Get list of SNPs from selected CpGs
snps_cojo <- unique(smr_noheidi_bonf[probeID %in% cpgs_cojo, ]$topSNP)
# 176 SNPs for conditional analysis
write.table(x = snps_cojo, file = "./temp_top/conditional_list_of_snps.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#############################################
#####################################
##
## 1st ROUND OF COJO: CONDITIONING ON TOP-SNP
##
#####################################
#############################################

system("gcta --bfile dubois_genotype --cojo-file ./PGC3_SCZ_wave3.extended.autosome.public.v3.freq.ma --cojo-cond ./temp_top/conditional_list_of_snps.txt --out ./temp_top/colocSCZ_genoDubois --cojo-wind 500 --cojo-slct")

# Upload COJO results
cojo_gwas <- fread("./temp_top/colocSCZ_genoDubois.cma.cojo")
table(cojo_gwas$pC < 5e-8)
# FALSE   TRUE 
# 167433     10 
cojo_gwas <- cojo_gwas[pC < 5e-8, ]

#############################################
#####################################
##
## 2nd ROUND OF COJO: CONDITIONING ON SECONDARY SIGNAL
##
#####################################
#############################################

### Get list of SNPs from selected CpGs
snps_cojo2 <- unique(cojo_gwas$SNP)
# 10 SNPs to use for conditioning the whole GWAS

write.table(x = snps_cojo2, file = "./temp_second/conditional_list_of_snps.txt", quote = F, sep = "\t", row.names = F, col.names = F)

system("gcta --bfile ./dubois_genotype --cojo-file ./PGC3_SCZ_wave3.extended.autosome.public.v3.freq.ma --cojo-cond ./temp_second/conditional_list_of_snps.txt --out ./temp_second/colocSCZ_genoDubois --cojo-wind 500 --cojo-slct")

# Upload results
cojo_gwas2 <- fread("./temp_second/colocSCZ_genoDubois.cma.cojo")
cojo_gwas2 <- cojo_gwas2[!is.na(pC), ]

# Upload complete GWAS 
gwas <- fread("./PGC3_SCZ_wave3.extended.autosome.public.v3.freq.ma")

# Check allele
table(cojo_gwas2[SNP %in% intersect(cojo_gwas2$SNP, gwas$SNP), ]$refA == gwas[SNP %in% intersect(cojo_gwas2$SNP, gwas$SNP),]$A1)
# TRUE 
# 173776 

# Add A2 info into COJO results
cojo_gwas2 <- merge(x = cojo_gwas2, y = gwas[,c("SNP","A1","A2")], by = "SNP")

# Write GWAS in .ma format
cojo_gwas2 <- cojo_gwas2[,c("SNP","A1","A2","freq","bC","bC_se","pC","n")]
colnames(cojo_gwas2) <- c("SNP","A1","A2","freq","b","se","p","n")
fwrite(x = cojo_gwas2, file = "./conditional_GWAS_SCZ.ma", quote = F, sep = "\t", row.names = F, col.names = T)
