tissue=$1
pheno=$2
for i in {1..22}
do
	Rscript ./FUSION.assoc_test.R --sumstats ../data/gwas/${pheno}.sumstats.gz --weights ../result/$tissue.pos --weights_dir ../result/ --ref_ld_chr ../data/LDREF/hg38/1000G.EUR.hg38. --chr $i --out ../assoc.result/$pheno.$tissue.ldrefhg38.${i} --force_model etwas --perm 10000
done
