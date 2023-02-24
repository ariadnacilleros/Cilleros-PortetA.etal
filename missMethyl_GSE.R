### Characterization of placental mQTLs with missMethyl
### Ariadna Cilleros-Portet

setwd("./")
library(missMethyl)
library(data.table)
library(ggplot2)

# TensorQTL results
nominal <- fread("nominalfastqtl_format_pval5e-8.txt", col.names = c("cpg","snp","dist","p","b"))

# Get the list of CpGs tested in mQTL mapping
bed <- fread("./nominalINMA_betas_0based.txt")
background_cpgs <- bed$ID
rm(bed)

# Gene Set Enrichment with GO
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
gst <- gometh(sig.cpg=unique(nominal$cpg), collection="GO", array.type = "EPIC", all.cpg = background_cpgs,
              plot.bias=TRUE)
table(gst$FDR < 0.05)
# FALSE  TRUE 
# 22334   340 
gst <- gst[gst$FDR < 0.05, ]
gst <- gst[order(gst$FDR, decreasing = T),]
gst$TERM <- factor(gst$TERM, levels=gst$TERM)

## Get enrichment barplot
gst_plot <- ggplot(gst[1:30, ], aes(x = TERM, y = N)) +
  geom_col(aes(fill=FDR)) + 
  scale_fill_gradient(low="red", high="blue")+
  guides(fill=guide_colorbar(title='FDR')) +
  labs(title="Gene Set Analysis with missMethyl GO gene sets",
       x ="GO terms", y = "Count")+
  coord_flip() + theme_bw()


# Gene Set Enrichment with KEGG
kegg <- gometh(sig.cpg=unique(nominal$cpg), collection="KEGG", array.type = "EPIC", all.cpg = background_cpgs,
              plot.bias=TRUE)
table(kegg$FDR < 0.05)
# FALSE  TRUE 
# 314    38 
kegg <- kegg[kegg$FDR < 0.05, ]
kegg <- kegg[order(kegg$FDR, decreasing = T),]
kegg$Description <- factor(kegg$Description, levels=kegg$Description)

## Get enrichment barplot
kegg_plot <- ggplot(kegg[1:30, ], aes(x = Description, y = N)) +
  geom_col(aes(fill=FDR)) + 
  scale_fill_gradient(low="red", high="blue")+
  guides(fill=guide_colorbar(title='FDR')) +
  labs(title="Gene Set Analysis with missMethyl KEGG gene sets",
       x ="KEGG terms", y = "Count")+
  coord_flip() + theme_bw()
