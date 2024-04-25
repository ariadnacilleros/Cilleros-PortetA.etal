## ReadMe A. Cilleros-Portet et al.2023
This repository contains the coded used for the A. Cilleros-Portet et al.2023, and the scripts are described as follows. 
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
