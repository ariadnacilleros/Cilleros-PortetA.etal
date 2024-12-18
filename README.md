## Potentially causal associations between placental DNA methylation and schizophrenia and other neuropsychiatric disorders

[![DOI](https://zenodo.org/badge/594034394.svg)](https://doi.org/10.5281/zenodo.14198404)

[**Ariadna Cilleros-Portet**](ariadna.cilleros@ehu.eus), Corina Lesseur, Sergi Marí, Marta Cosin-Tomas, Manuel Lozano, Amaia Irizar, Amber Burt, Iraia García-Santisteban, Diego Garrido Martín, Geòrgia Escaramís, Alba Hernangomez-Laderas, Raquel Soler-Blasco, Charles E Breeze, Bárbara P Gonzalez-Garcia, Loreto Santa-Marina, Jia Chen, Sabrina Llop, Mariana F Fernández, Martine Vrijhed, Jesús Ibarluzea, Mònica Guxens, Carmen Marsit, Mariona Bustamante, Jose Ramon Bilbao, [**Nora Fernandez-Jimenez**](nora.fernandez@ehu.eus)

Department of Genetics, Physical Anthropology and Animal Physiology, Biocruces-Bizkaia Health Research Institute and University of the Basque Country (UPV/EHU), Leioa, Spain.

Increasing evidence supports the role of placenta in neurodevelopment and potentially, in the later onset of neuropsychiatric disorders. Recently, methylation quantitative trait loci (mQTL) and interaction QTL (iQTL) maps have proven useful to understand SNP-genome wide association study (GWAS) relationships, otherwise missed by conventional expression QTLs. In this context, we propose that part of the genetic predisposition to complex neuropsychiatric disorders acts through placental DNA methylation (DNAm). We constructed the first public placental cis-mQTL database including nearly eight million mQTLs calculated in 368 fetal placenta DNA samples from the INMA project, ran cell type- and gestational age-imQTL models and combined those data with the summary statistics of the largest GWAS on 10 neuropsychiatric disorders using Summary-based Mendelian Randomization (SMR) and colocalization. Finally, we evaluated the influence of the DNAm sites identified on placental gene expression in the RICHS cohort. We found that placental cis-mQTLs are highly enriched in placenta-specific active chromatin regions, and useful to map the etiology of neuropsychiatric disorders at prenatal stages. Specifically, part of the genetic burden for schizophrenia, bipolar disorder and major depressive disorder confers risk through placental DNAm. The potential causality of several of the observed associations is reinforced by secondary association signals identified in conditional analyses, regional pleiotropic methylation signals associated to the same disorder, and cell type-imQTLs, additionally associated to the expression levels of relevant immune genes in placenta. In conclusion, the genetic risk of several neuropsychiatric disorders could operate, at least in part, through DNAm and associated gene expression in placenta.

Placental mQTLs browser: ()[]

## ReadMe 
This repository contains the coded used for the A. Cilleros-Portet et al. 2023, and the scripts are described as follows. 
1. Map of cis-mQTLs with TensorQTL nominal approach ([tensorqtl_nominal.py](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/tensorqtl_nominal.py)). 
2. Map of cis-mQTLs with TensorQTL permuted and conditional approach ([tensorqtl_permuted_conditional.py](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/tensorqtl_permuted_conditional.py)).
3. Map of interacting cis-mQTLs with TensorQTL interaction approach ([tensorqtl_interaction.py](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/tensorqtl_interaction.py)).
4. Get the final mQTL database, only applicable for nominal and interaction analysis ([get_final_database.sh](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/get_final_database.sh)).
5. Gene set enrichment analysis with Disease Onthology database ([GSEA_analysis.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/GSEA_analysis.R)).
6. Over-representation analysis with Gene Onthology and Kyoto Encyclopedia of Genes and Genomes ([missMethyl_GSE.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/missMethyl_GSE.R)).
7. Illumina annotation enrichment chi-square tests ([enrichment.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/enrichment.R)).
8. Select CpG list for eFORGE analysis ([CpGlist_eFORGE.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/CpGlist_eFORGE.R)).
9. Formatting GWAS summary-statistics for SMR ([checkalleles_maformat.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/checkalleles_maformat.R)).
10. Formatting mQTL database for SMR ([getFastQTL_format.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/getFastQTL_format.R)).
11. Run mutli-SNP based SMR test ([runSMR.sh](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/runSMR.sh)).
12. Get the CpG list overlapping GWAS loci for colocalization analysis ([coloc_step1_overlap_CpGs_gr.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/coloc_step1_overlap_CpGs_gr.R))
13. Get the BED file for colocalization analysis ([coloc_step2_get_bed_file.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/coloc_step2_get_bed_file.R))
14. Run TensorQTL for colocalization analysis ([coloc_step3_run_tensorqtl.py](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/coloc_step3_run_tensorqtl.py))
15. Format TensorQTL results for colocalization analysis ([coloc_step4_concat_tensor_results.sh](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/coloc_step4_concat_tensor_results.sh))
16. Run colocalization test ([coloc_step5_run_coloc.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/coloc_step5_run_coloc.R))
17. Classify cell-type interacting-mQTLs according to [Kim-Hellmuth et al.2020](https://www.science.org/doi/10.1126/science.aaz8528) interpretation ([classify_imQTLs.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/classify_imQTLs.R)).
18. Condition analysis of the GWAS summary-statistics ([conditional_GWAS.r](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/conditional_GWAS.r)).
19. Condition analysis of the mQTLs summary-statistics ([conditional_mQTLs.r](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/conditional_mQTLs.r)). 
20. Run multi-SNP based SMR test using the conditional GWAS and mQTL summary-statistics ([conditional_SMR.R](https://github.com/ariadnacilleros/Cilleros-PortetA.etal/blob/main/conditional_SMR.R)).
