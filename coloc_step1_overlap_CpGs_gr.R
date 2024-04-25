#### Overlap between CpGs permuted FDR < 0.05 and coloc regions 
#### Ariadna Cilleros-Portet

setwd("./")
library(data.table)
library(GenomicRanges)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(stringr)

# Read permuted database 
perm <- fread("./cis_tensorQTL_maf05_hwe05_PC5_sex_planet_18mPCresidual_RNT_PERMUT30102023_FDR.txt")
head(perm)
length(perm$phenotype_id)
# 214,830 CpGs

## Add CpG position and chromosome
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
perm <- merge(x = perm, y = annot[,c("Name","chr","pos")], by.x = "phenotype_id", by.y = "Name")
head(perm)
perm$chr <- str_remove(perm$chr, pattern = 'chr')

# Read SCZ regions
regions <- data.table(readxl::read_xls("./PGC3_SCZ_wave3.extended.287regions.xls"))
regions$`merge-LEFT` <- as.numeric(regions$`merge-LEFT`)
regions$`merge-RIGHT` <- as.numeric(regions$`merge-RIGHT`)
regions$Chromosome <- as.numeric(regions$Chromosome)
head(regions)

# Add 500kb per side 
regions$`merge-LEFT-500` <- regions$`merge-LEFT`-500000
regions$`merge-RIGHT-500` <- regions$`merge-RIGHT`+500000

# Create GenomicRanges object 
regions_gr <- makeGRangesFromDataFrame(df = regions,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="Chromosome",
                         start.field="merge-LEFT-500",
                         end.field="merge-RIGHT-500")

perm_gr <- makeGRangesFromDataFrame(df = perm,
                         keep.extra.columns=T,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="pos", end.field = "pos")

# Finding overlap 
findOverlaps(perm_gr, regions_gr)
perm_overlap <- perm_gr[queryHits(findOverlaps(perm_gr, regions_gr, type = 'any')),]
# 38,412 from 214,830 CpGs
write.table(x = perm_overlap$phenotype_id, file = "./list_of_CpG_coloc.txt", quote = F, sep = "\t", row.names = F, col.names = T)
