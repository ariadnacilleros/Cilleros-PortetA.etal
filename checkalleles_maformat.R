### Prepare GWAS for SMR (.ma format)
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)
library(stringr)

# Load GWAS full summary-statistics file)
gwas <- fread("./PGC3_SCZ_wave3.extended.autosome.public.v3.vcf.tsv")

# Load genotype information date
bim <- fread("./INMA.bim")

# Create chr:pos id
gwas$snp <- paste(gwas$CHROM, gwas$POS, sep=":")

# Create N column
gwas$n <- gwas$NCAS + gwas$NCON

# Filter GWAS and genotype files
bimsnps <- intersect(bim[, V2], gwas$snp)
gwas <- gwas[match(bimsnps, snp), .(snp, A1, A2, FCON, BETA, SE, PVAL, n)]
bim <- bim[V2%in%bimsnps, ]

# Adjust effect allele between GWAS and genotype files
gwas[A1 == bim$V6, c('A1', 'A2', 'FCON', 'BETA') := .(A2, A1, 1-FCON, -BETA)]

# Change colnames according .ma format
colnames(gwas) <- c("SNP","A1","A2","freq","b","se","p","n")

# Write output file
fwrite(x = gwas, file = "./PGC3_SCZ_wave3.extended.autosome.public.v3.freq.ma", sep = "\t", 
       quote = F, col.names = T, row.names = F, na = "NA")
