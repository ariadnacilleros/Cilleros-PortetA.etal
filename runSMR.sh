### Get file set and run SMR
### Ariadna Cilleros-Portet

# Get SMR mQTL summary
smr --eqtl-summary ./nominal_db_filterpval_fastqtl.txt --fastqtl-nominal-format --make-besd --out nominal_db

# Run multi-SMR 
smr --bfile ./inma --gwas-summary ./PGC3_SCZ_wave3.extended.autosome.public.v3.freq.ma --beqtl-summary ./nominal_db --out nominal_inma_outcomeSCZ_multi --smr-multi --set-wind 500 --ld-multi-snp 0.9

