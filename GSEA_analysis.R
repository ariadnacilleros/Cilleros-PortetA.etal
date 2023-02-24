### Characterization of placental mQTLs with GSEA
### Ariadna Cilleros-Portet

setwd("./")
library(data.table)
library(clusterProfiler)
library(ggplot2)
library(DOSE)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(methylGSA)
library(tibble)
library(dplyr)

# TensorQTL results
nominal <- fread("nominalfastqtl_format_pval5e-8.txt", col.names = c("cpg","sno","dist","p","b"))

### GSEA Disease
# Get genes 
annot <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annot <- annot[nominal$cpg, ]
genes <- annot[annot$UCSC_RefGene_Name !="", ]
out <- do.call(rbind, Map(cbind, strsplit(genes$UCSC_RefGene_Name, ";"), genes$Name))
out <- as_tibble(out)
out2 <- out %>% distinct()

# Get gene ID from symbol to entrezID
genes_gsea <- merge(x = out2, y = nominal[,c("cpg", "p")], by.x="V2", by.y="cpg")
genes_entrez <- bitr(genes_gsea$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genes_gsea_entrez <- merge(x = genes_entrez, y = genes_gsea, by.x= "SYMBOL", by.y="V1")

# Get named vector with p-value and CpGs id
geneList <- as.list(genes_gsea_entrez$p)
names(geneList) <- c(genes_gsea_entrez$ENTREZID)

geneList <- setNames(genes_gsea_entrez$p, genes_gsea_entrez$ENTREZID)

## Disease onthology
DO <- gseDO(geneList = sort(geneList, decreasing = T),
            minGSSize     = 120,
            #pvalueCutoff  = 5e-8,
            #pAdjustMethod = "bonferroni",
            verbose       = FALSE)
do <- setReadable(DO, 'org.Hs.eg.db')
do <- as.data.frame(do)

### Get enrichment dotplot
gp <- enrichplot::dotplot(object = DO, showCategory=30) + ggtitle("GSEA Disease onthology")
ggsave(filename = "./methylgsea_dotplot.png", device = "png", plot = gp, dpi = 320, width = 10, height = 10)

