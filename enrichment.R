### Characterization of placental mQTLs with annotation
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)

# TensorQTL results
nominal <- fread("nominalfastqtl_format_pval5e-8.txt", col.names = c("cpg","sno","dist","p","b"))

# Get annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Get annotation df for sign and non-sign CpGs
annot_sign <- annot[annot$Name %in% unique(nominal$cpg), ]
annot_notsign <- annot[annot$Name %in% setdiff(unique(annot$Name), nominal$cpg), ]

## Enrichment of Exon Bond
length(annot_sign[grep("ExonBnd",annot_sign$UCSC_RefGene_Group),]$Name)
length(annot_notsign[grep("ExonBnd",annot_notsign$UCSC_RefGene_Group),]$Name)

#create a dataframe for chi-square test, with two columns: significatives, and Exon Bond overlap
sign_exon <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                            exon = c(rep("yes", 747), rep("no", 109974), 
                                         rep("yes", 6499), rep("no",748639)))
chisq_sign_exon <- chisq.test(sign_exon$sign, sign_exon$exon)
round(chisq_sign_exon$residuals, 3)

## Enrichment of 3'UTR
length(annot_sign[grep("3'UTR",annot_sign$UCSC_RefGene_Group),]$Name)
length(annot_notsign[grep("3'UTR",annot_notsign$UCSC_RefGene_Group),]$Name)

#create a dataframe for chi-square test, with two columns: significatives, and 3' UTR overlap
sign_3utr <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                        utr = c(rep("yes", 2918), rep("no", 107803), 
                                 rep("yes", 21538), rep("no",733600)))
chisq_sign_3utr <- chisq.test(sign_3utr$sign, sign_3utr$utr)
round(chisq_sign_3utr$residuals, 3)


## Enrichment of promoter regions
length(annot_sign[grep("TSS1500|TSS200|5'UTR|1stExon",annot_sign$UCSC_RefGene_Group),]$Name)
length(annot_notsign[grep("TSS1500|TSS200|5'UTR|1stExon",annot_notsign$UCSC_RefGene_Group),]$Name)

#create a dataframe for chi-square test, with two columns: significatives, and promoter overlap
sign_promoter <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                            promoter = c(rep("yes", 29966), rep("no", 80755), 
                                         rep("yes", 262932), rep("no",492206)))
chisq_sign_promoter <- chisq.test(sign_promoter$sign, sign_promoter$promoter)
round(chisq_sign_promoter$residuals, 3)

## Enrichment of body regions
length(annot_sign[grep("Body",annot_sign$UCSC_RefGene_Group),]$Name)
length(annot_notsign[grep("Body",annot_notsign$UCSC_RefGene_Group),]$Name)

#create a dataframe for chi-square test, with two columns: significatives, and body
sign_body <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                            body = c(rep("yes", 46991), rep("no", 63730), 
                                         rep("yes", 313728), rep("no",441410)))
chisq_sign_body <- chisq.test(sign_body$sign, sign_body$body)
round(chisq_sign_body$residuals, 3)

# Enrichment of islands / shelf / shore / open sea 
table(annot_sign$Relation_to_Island)
table(annot_notsign$Relation_to_Island)

#create a dataframe for chi-square test, with two columns: significatives, and Island
sign_island <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                            island = c(rep("yes", 9509), rep("no", 101212), 
                                         rep("yes", 151932), rep("no",603206)))
chisq_sign_island <- chisq.test(sign_island$sign, sign_island$island)
round(chisq_sign_island$residuals, 3)



#create a dataframe for chi-square test, with two columns: significatives, and shelf 
sign_shelf <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                          shelf = c(rep("yes", 8295), rep("no", 102426), 
                                     rep("yes", 53396), rep("no",701742)))
chisq_sign_shelf <- chisq.test(sign_shelf$sign, sign_shelf$shelf)
round(chisq_sign_shelf$residuals, 3)
# enrichment

#create a dataframe for chi-square test, with two columns: significatives, and shore 
sign_shore <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                         shore = c(rep("yes", 21266), rep("no", 89455), 
                                   rep("yes", 133280), rep("no",621858)))
chisq_sign_shore <- chisq.test(sign_shore$sign, sign_shore$shore)
round(chisq_sign_shore$residuals, 3)
# enrichment

#create a dataframe for chi-square test, with two columns: significatives, and open sea
sign_opensea <- data.frame(sign = c(rep("yes", 110721),rep("no", 755138)), 
                         opensea = c(rep("yes", 71651), rep("no", 39070), 
                                   rep("yes", 416530), rep("no",338608)))
chisq_sign_opensea <- chisq.test(sign_opensea$sign, sign_opensea$opensea)
round(chisq_sign_opensea$residuals, 3)
# enrichment

