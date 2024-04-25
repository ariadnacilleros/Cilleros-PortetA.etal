### Get BED file for coloc SCZ CpG list 
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

## Load CpG list
cpgs <- fread("./list_of_CpG_coloc.txt")
head(cpgs)
#38,412 CpGs
length(unique(cpgs$x))
# 35,996 unique CpGs

## Load BED file 
epic <- fread("./nominalINMA_ranknormal_0based.txt")
head(epic)
dim(epic)

## Filter BED file by CpG list
table(cpgs$x %in% epic$ID) 
# TRUE 
# 38412 
epic_filt <- epic[ID %in% cpgs$x, ]
dim(epic_filt)
# 35,996 CpGs
head(epic_filt)

## Write BED file
fwrite(x = epic_filt, file = "./nominalINMA_ranknormal_0based_colocSCZ.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, )

## Sort BED file
system("(head -n1 nominalINMA_ranknormal_0based_colocSCZ.txt && sort -k1,1V -k2,2n -k3,3n <(tail -n+2 nominalINMA_ranknormal_0based_colocSCZ.txt)) > nominalINMA_ranknormal_0based_colocSCZ_sorted.bed")

## Zip BED file
system("bgzip nominalINMA_ranknormal_0based_colocSCZ_sorted.bed")
