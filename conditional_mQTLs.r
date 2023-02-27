### Conditional signals on mQTLs
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# Obtain mQTLs text file (FastQTL format) with the significant conditional signals (and its statistics)
# Load coloc results
coloc <- as.data.table(readRDS(file = "./output_coloc.RDS"))

# Subset coloc results by PPA4 > 0.8 (ref. Leyden et al.2022)
coloc <- coloc[coloc$PP.H4.abf > 0.8, ]
#2,177 coloc hits; 1603 CpGs, 188 regions

# Load SMR+HEIDI results
smr <- fread("./nominal_inma_outcomeSCZ_multi.msmr")

# Get COJO hits
smr$bonf <- p.adjust(smr$p_SMR_multi, "bonferroni")
smr_noheidi_bonf <- smr[bonf < 0.05, ]
smr_noheidi_bonf <- smr_noheidi_bonf[p_HEIDI < 0.05, ]
# 469 SMR hits with heterogeneity (p-HEIDI < 0.05)

# Intersect SMR hits with heterogeneity with coloc to perform COJO
cpgs_cojo <- intersect(smr_noheidi_bonf$probeID, coloc$trait2)
# 243 CpGs con evidencia de heterogeniedad, significativas en SMR y con evidencia de colocalización 

# Get mQTLs data in .ma format to perform cojo in mQTLs database 
mqtls <- fread("./nominal_db_filterpval_fastqtl.txt", col.names = c("phenotype_id","variant_id","tss_distance","af","ma_samples","ma_count","pval_nominal","slope","se"))
mqtls_subset <- mqtls[phenotype_id %in% cpgs_cojo, ]
rm(mqtls)

# Load alleles from bim file
bim <- fread("./inma.bim")
colnames(bim)[5:6] <- c("A1","A2")

# Add alleles information to mQTLs data
mqtls_subset <- merge(x = mqtls_subset, y = bim[,c("V2","A1","A2")], by.x="variant_id", by.y="V2")

#############################################
#####################################
##
## 1st ROUND OF COJO: CONDITIONING ON TOP-SNP
##
#####################################
#############################################

# For loop per CpG 
for (cpg in cpgs_cojo){
  
  # Write the text file with the top SNP (by SMR) for the corresponding CpG
  write.table(x = smr_noheidi_bonf[probeID == cpg, ]$topSNP, file = "./temp_conditional_list_of_snps.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Subset nominal mQTLs database by the corresponding CpG
  mqtls_cpg <- mqtls_subset[phenotype_id == cpg, ]
  mqtls_cpg_ma <- mqtls_cpg[,c("variant_id","A1","A2","af","slope","se","pval_nominal")]
  mqtls_cpg_ma$N <- 368
  colnames(mqtls_cpg_ma) <- c("SNP","A1","A2","freq","b","se","p","N")
  
  # Write in .ma format the subseted nominal mQTLs database 
  fwrite(x = mqtls_cpg_ma, file = "./temp_nominal_inma.ma", quote = F, row.names = F, sep = "\t")
  
  # Run GCTA cojo analysis for the corresponding CpG
  system("gcta --bfile ./dubois_genotype --cojo-file ./temp_nominal_inma.ma --cojo-cond ./temp_conditional_list_of_snps.txt --out ./temp_mQTLs_colocnominal_inma_genoDubois --cojo-wind 500 --cojo-slct")
  
  # Upload results & write a text file with the CpG id
  temp_results <- fread("./temp_mQTLs_colocnominal_inma_genoDubois.cma.cojo", header = T)
  write.table(x = rep(smr_noheidi_bonf[probeID == cpg, ]$probeID,nrow(temp_results)), file = "./temp_cpg.txt", quote = F, row.names = F, col.names = T, sep = "\t")
  
  # Paste CpG id with conditional results
  system("paste temp_cpg.txt temp_mQTLs_colocnominal_inma_genoDubois.cma.cojo >> all_results.txt") 
}

# Upload results
allresults <- fread("./all_results.txt")
allresults <- allresults[x != "x", ]
length(unique(allresults$x))
#240 CpGs

# Hits with evidence of secundary signals 
allresults$pC <- as.numeric(allresults$pC)
table(allresults$pC < 5e-8) 
# FALSE  TRUE 
# 47890   332
allresults <- allresults[pC < 5e-8, ]

#############################################
#####################################
##
## 2nd ROUND OF COJO: CONDITIONING ON SECONDARY SIGNAL
##
#####################################
#############################################

# Get mQTLs data in .ma format to perform cojo in mQTLs database 
mqtls <- fread("./nominal_db_filterpval.txt", col.names = c("phenotype_id","variant_id","tss_distance","af","ma_samples","ma_count","pval_nominal","slope","se"))
mqtls_subset <- mqtls[phenotype_id %in% unique(allresults$x), ]
rm(mqtls)

# Load alleles from bim file
bim <- fread("./inma.bim")
colnames(bim)[5:6] <- c("A1","A2")

# Add alleles information to mQTLs data
mqtls_subset <- merge(x = mqtls_subset, y = bim[,c("V2","A1","A2")], by.x="variant_id", by.y="V2")

# For loop per CpG 
for (cpg in unique(allresults$x)){
  
  # Write the text file with the conditional SNP for the corresponding CpG
  write.table(x = allresults[x == cpg, ]$SNP, file = "./temp_second/temp_conditional_list_of_snps.txt", quote = F, row.names = F, col.names = F, sep = "\t")
  
  # Subset nominal mQTLs database by the corresponding CpG
  mqtls_cpg <- mqtls_subset[phenotype_id == cpg, ]
  mqtls_cpg_ma <- mqtls_cpg[,c("variant_id","A1","A2","af","slope","se","pval_nominal")]
  mqtls_cpg_ma$N <- 368
  colnames(mqtls_cpg_ma) <- c("SNP","A1","A2","freq","b","se","p","N")
  
  # Write in .ma format the subseted nominal mQTLs database 
  fwrite(x = mqtls_cpg_ma, file = "./temp_second/temp_nominal_inma.ma", quote = F, row.names = F, sep = "\t")
  
  # Run GCTA cojo analysis for the corresponding CpG
  system("gcta --bfile ./dubois_genotype --cojo-file ./temp_second/temp_nominal_inma.ma --cojo-cond ./temp_second/temp_conditional_list_of_snps.txt --out ./temp_second/temp_mQTLs_colocnominal_inma_genoDubois --cojo-wind 500 --cojo-slct")
  
  # Upload results & write a text file with the CpG id
  temp_results <- fread("./temp_mQTLs_colocnominal_inma_genoDubois.cma.cojo", header = T)
  write.table(x = rep(cpg,nrow(temp_results)), file = "./temp_second/temp_cpg.txt", quote = F, row.names = F, col.names = T, sep = "\t")
  
  # Paste CpG id with conditional results
  system("paste temp_second/temp_cpg.txt temp_second/temp_mQTLs_colocnominal_inma_genoDubois.cma.cojo >> temp_second/all_results.txt") 
}

# Upload results from COJO secondary signal
allresults2 <- fread("./temp_second/all_results.txt")
allresults2 <- allresults2[x != "x", ]
allresults2 <- allresults2[!is.na(pC), ] #labelled as "NA" if the multivariate correlation between the SNP in question and all the covariate SNPs is > 0.9
allresults2$mqtl <- paste(allresults2$x, allresults2$SNP, sep="-")

# Upload nominal mQTLs db 
nominal <- fread("./nominal_db_filterpval.txt", col.names = c("phenotype_id","variant_id","dist","pval","b"))
nominal$mqtl <- paste(nominal$phenotype_id, nominal$variant_id, sep="-")

# Add distance CpG-SNP 
allresults2 <- merge(x = allresults2, y = nominal[,c("dist","mqtl")], by="mqtl")

# Create FastQTL format for conditional mQTLs database 
allresults2_fastqtl <- allresults2[,c("x","SNP","dist","pC","bC")]
fwrite(x = allresults2_fastqtl, file = "./conditional_mQTLs_db.txt", quote = F, sep = "\t", col.names = F, row.names = F)

