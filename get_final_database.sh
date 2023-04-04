## Get the final nominal/interaction mQTL database from the text files 
## Ariadna Cilleros-Portet

#!/bin/bash
rm nominal_db_filterpval.txt
touch nominal_db_filterpval.txt

# Filter each chromosome txt by p-value < 5e-8, and write a unique file with the mQTLs from all chromosomes
for i in {1..22}
	do
		awk '{ if ($7 <= 5e-8) print $0 }' ./INMA_nominal.${i}.txt >> nominal_db_filterpval.txt 
	done
