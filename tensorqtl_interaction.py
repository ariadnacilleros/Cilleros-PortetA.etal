## Run TensorQTL for the interaction database
## Ariadna Cilleros-Portet

# Load packages
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# Define paths to data
plink_prefix_path = './plink_prefix_path'
expression_bed = './expression_bed.bed.gz'
covariates_file = './covariates_file.txt'
prefix = './INMA_interaction'

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

# Get interaction term
import numpy as np
ga_pd = pd.read_csv('./interact.txt', sep='\t', index_col=0, header=0)
ga_array = np.array(ga_pd['Syncytiotrophoblast'])
interaction_s = pd.Series(ga_array,index=ga_pd.index)
interaction_s.sort_index(inplace=True)

# Run TensorQTL interaction approach
cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix,
                covariates_df=covariates_df,
                interaction_s=interaction_s, 
				window=500000, 
				output_dir='./', write_top=True, write_stats=True)

# Write .parquet to .txt
for x in range(1,23):
	pairs_df = pd.read_parquet(f'{prefix}.cis_qtl_pairs.{x}.parquet')
	pairs_df.to_csv(f'./interaction_database.{x}.txt', header=True, index=False, sep='\\t')
