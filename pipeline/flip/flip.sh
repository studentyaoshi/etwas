#  The pipeline is used to align all alleles on the forward strand.

## Download all the genome files in ~/etwas/pipeline/flip/chr{1..22}.fa

## Extract SNPs without ambiguity
python checkstring.py data.bim > withchain.bim  
awk -F '\t' {'print$2'} withchain.bim> snp.list

## Generate the genotype files with only SNPs without ambiguity
plink --bfile data --extract snp.list --make-bed --out data2
python new_bim.py withchain.bim data2.bim > data3.bim
mv data3.bim data2.bim

## Extract SNPs on the reverse strand.
awk -F '\t' '$7=="-"{print$2}' withchain.bim > minus.snp

## Flip the strand of allele on the reverse strand.
plink --bfile data2 --flip minus.snp --make-bed --out data_flip
