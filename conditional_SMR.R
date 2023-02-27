### Run SMR analysis on conditional signals
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# Get .besd, .epi, .esi files 
system("smr --eqtl-summary ./conditional_mQTLs_db.txt --fastqtl-nominal-format --make-besd --out conditional_mQTLs_db")

# Perform SMR with both conditional signals  
system("smr --bfile ./inma --gwas-summary ./conditional_GWAS_SCZ.ma --beqtl-summary ./conditional_mQTLs_db --out ./bothconditionalsignals --smr-multi --set-wind 500 --ld-multi-snp 0.9")
condBoth <- fread("./bothconditionalsignals.msmr")
condBoth$bonf <- p.adjust(p = condBoth$p_SMR_multi, method = "bonferroni")
condBoth <- condBoth[bonf < 0.05, ]
dim(condBoth[p_HEIDI > 0.05 | is.na(p_HEIDI), ])
