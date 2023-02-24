## Run TensorQTL for the nominal database
## Ariadna Cilleros-Portet

# Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print('PyTorch {}'.format(torch.__version__))
print('Pandas {}'.format(pd.__version__))

# Define paths to data
plink_prefix_path = './wg_filt_maf05_hwe05_multiallelic_chrpos'
expression_bed = './nominalINMA_betas_0based.bed.gz'
covariates_file = './maf05_hwe05_COV_sex_PC5_planet.txt'

# Load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Set column order
phenotype_df = phenotype_df.reindex(sorted(phenotype_df.columns), axis=1)
genotype_df = genotype_df.reindex(sorted(genotype_df.columns), axis=1)
covariates_df = covariates_df.sort_index()

# Run TensorQTL permuted approach
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df, seed=123456789)

# Run TensorQTL conditional approach 
indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                               phenotype_df, phenotype_pos_df, covariates_df)