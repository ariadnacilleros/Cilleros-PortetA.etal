# Run TensorQTL
# Ariadna Cilleros-Portet
# python 3.6

import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# Adjust the number of CPU cores here
num_threads = 8 
torch.set_num_threads(num_threads)

# Define paths to data
plink_prefix_path = './wg_filt_maf05_hwe05_multiallelic_chrpos'
expression_bed = './nominalINMA_ranknormal_0based_colocSCZ_sorted.bed.gz'
covariates_file = './maf05_hwe05_RNT_sex_5PC_planet_18mPCresidual.txt'
prefix = './parquet/INMA_nominal_260324'

# Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Sort columns (samples)
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

# Run TensorQTL
cis.map_nominal(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix, covariates_df=covariates_df, 
                window=500000)
                
# Transform .parquet to text file
for x in range(1,23):
    pairs_df = pd.read_parquet(f'{prefix}.cis_qtl_pairs.{x}.parquet')
    print(x)
    pairs_df.to_csv(f'./chr_txt/INMA_nominal_260324.{x}.txt', header=True, index=False, sep='\t') 
