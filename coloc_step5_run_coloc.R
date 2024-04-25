#### Run coloc analysis with SCZ for nominal model
#### Ariadna Cilleros-Portet

setwd("./")
library(data.table)
library(coloc)
library(stringr)

# Load mQTLs database
nominal <- fread("chr_txt/INMA_nominal_260324.all_chr.txt")
head(nominal)
# 57,081,448 mQTLs
nominal <- nominal[, maf := ifelse(1-af < af, 1-af, af)]

# Load CpG annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Add CpG position in mQTLs database
nominal <- merge(x = nominal, y = annot[,c("chr","pos","Name")], by.x = "phenotype_id", by.y = "Name")
nominal$chr <- as.numeric(str_replace(string = nominal$chr, pattern = "chr", replacement = ""))
nominal$pos <- as.numeric(nominal$pos)

# Load regions
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
          dataset1<-list(beta=regions_mqtls$slope,varbeta=(regions_mqtls$slope_se)^2,type = "quant",snp = regions_mqtls$variant_id, MAF = regions_mqtls$maf, N=368) ### mQTL results
          
          # Run coloc
          my.res <- coloc.abf(dataset1, dataset2)
          output<-rbind(output, c("region"=paste(chr_regions$`merge-RIGHT`[region],"-",chr_regions$`merge-LEFT`[region]), "trait2"=probe, my.res$summary))
        }
      }
    }
  }
}

saveRDS(object = output, file = "./output_coloc.RDS")
output <- as.data.table(readRDS(file = "./output_coloc.RDS"))

output$`PP.H3.abf+PP.H4.abf` <- as.numeric(output$PP.H3.abf) + as.numeric(output$PP.H4.abf)
dim(output[as.numeric(output$`PP.H3.abf+PP.H4.abf`) > 0.9, ])
# 21,071 hits

# Add lead SNP from genomic region (GWAS)
head(output)
head(regions)
regions$Region_coloc <- paste(regions$`merge-RIGHT`, regions$`merge-LEFT`, sep=" - ")
output2 <- merge(x = output, y = regions[,c("top-index","Region_coloc")], by.y = "Region_coloc", by.x="region")

# Add bp and chromosome
output2 <- merge(x = output2, y = annot[,c("Name", "chr", "pos")], by.x = "trait2", by.y = "Name")

# Order region coordinates 
library(tidyr)
output2 <- separate(data = as.tibble(output2), col = region, sep = " - ", into = c("end", "start")) 
output2$region <- paste(output2$start, output2$end, sep ="-")

# Filter and order columns
colnames(output2)
colnames(output2)[10] <- "PP.H3.abf+PP.H4.abf"
output2 <- output2[,c("chr", "region", "top.index", "nsnps", "trait2", "pos", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf","PP.H3.abf+PP.H4.abf")]
library(dplyr)
output2 <- as.data.table(arrange(output2, chr,region))
write.csv(x = output2[as.numeric(output2$`PP.H3.abf+PP.H4.abf`) > 0.9, ], file = "./sdata_coloc_scz.csv", quote = F, row.names = F)
openxlsx::write.xlsx(x = output2[as.numeric(output2$`PP.H3.abf+PP.H4.abf`) > 0.9, ], file = "./sdata_coloc_scz.xlsx", quote = F)
