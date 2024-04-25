ls chr_txt/* > input.txt
awk 'NR==1 || FNR > 1' $(< input.txt ) >> chr_txt/INMA_nominal_260324.all_chr.txt