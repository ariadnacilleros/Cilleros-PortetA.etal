#### Run coloc analysis with GWAS and nominal mQTLs db
#### Ariadna Cilleros-Portet

setwd("./")
library(data.table)
library(coloc)
library(stringr)

# Load mQTLs database
nominal <- fread("./nominal_db_filterpval.txt", col.names = c("phenotype_id","variant_id","tss_distance","af","ma_samples","ma_count","pval_nominal","slope","se"))
nominal <- nominal[, maf := ifelse(1-af < af, 1-af, af)]

# Load CpG annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Add CpG position in mQTLs database
nominal <- merge(x = nominal, y = annot[,c("chr","pos","Name")], by.x = "phenotype_id", by.y = "Name")
nominal$chr <- as.numeric(str_replace(string = nominal$chr, pattern = "chr", replacement = ""))
nominal$pos <- as.numeric(nominal$pos)

# Load region
regions <- data.table(readxl::read_xls("./PGC3_SCZ_wave3.extended.287regions.xls"))
regions$`merge-LEFT` <- as.numeric(regions$`merge-LEFT`)
regions$`merge-RIGHT` <- as.numeric(regions$`merge-RIGHT`)
regions$Chromosome <- as.numeric(regions$Chromosome)

# Load GWAS data
gwas <- fread("./PGC3_SCZ_wave3.extended.autosome.public.v3.ma")
gwas$CHR <- as.numeric(sub(x = gwas$SNP, pattern = ":.*", replacement = ""))
temp  <- str_split_fixed(string = gwas$SNP, pattern = ":", n = 2)
gwas$bp <- as.numeric(temp[,2])
rm(temp)
gwas <- gwas[, maf := ifelse(1-freq < freq, 1-freq, freq)]

# Coloc loop
output <- NULL
for (chr in as.numeric(unique(nominal$chr))){
  chr_mqtls <- nominal[nominal$chr == chr, ]
  chr_regions <- regions[regions$Chromosome == chr, ]
  
  for (region in 1:nrow(chr_regions)){
    regions_gwas<-gwas[which(gwas$CHR == chr & gwas$bp <= chr_regions$`merge-RIGHT`[region]+500000 & gwas$bp >= chr_regions$`merge-LEFT`[region]-500000),]
    if (nrow(regions_gwas) > 0){
      dataset2<-list(beta= regions_gwas$b, varbeta= (regions_gwas$se^2), type = "cc", snp = regions_gwas$SNP, MAF = regions_gwas$maf, N=regions_gwas$n) ### for each region SCZ GWAS results stay the same
      probes<-as.character(nominal[which(nominal$chr == chr & nominal$pos <= chr_regions$`merge-RIGHT`[region]+500000 & nominal$pos >= chr_regions$`merge-LEFT`[region]-500000),"phenotype_id"])
      
      for (probe in unique(probes)){
        regions_mqtls <- chr_mqtls[chr_mqtls$phenotype_id == probe, ] 
        if (nrow(regions_mqtls) > 1 & isTRUE(length(intersect(dataset2$snp, regions_mqtls$variant_id))>0)){
          dataset1<-list(beta=regions_mqtls$slope,varbeta=(regions_mqtls$se)^2,type = "quant",snp = regions_mqtls$variant_id, MAF = regions_mqtls$maf, N=368) ### mQTL results
          
          # Run coloc
          my.res <- coloc.abf(dataset1, dataset2)
          output<-rbind(output, c("region"=paste(chr_regions$`merge-RIGHT`[region],"-",chr_regions$`merge-LEFT`[region]), "trait2"=probe, my.res$summary))
        }
      }
    }
  }
}

# Get colocalization results (PPA4 > 0.8)
output[PP.H4.abf > 0.8, ]