tissue=$1
cat ../data/expression/gtex.gene.final.nomhc.anno|while read line
do gene=`echo $line|awk '{print$1}'`
	chr=`echo $line|awk '{print$2}'`
	position1=`echo $line|awk '{print$3}'`
	position2=`echo $line|awk '{print$4}'`
	let start=$position1-1000000
	let end=$position2+1000000
	if [ $start -lt 0 ]; then start=0; fi
	plink --bfile ../data/genotype/GTEx.chr$chr --maf 0.05 --chr $chr --from-bp $start --to-bp $end --make-bed --out ../tem/$gene
	gcta64 --bfile ../tem/$gene --make-grm --out ../tem/${gene} --thread-num 5
	mpheno_number=`head -1 ../data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.${tissue}.final|sed 's/\t/\n/g'|sed '1,2d' |grep -n -w ${gene}|awk -F ':' '{print$1}'`
	gcta64 --reml --grm ../tem/${gene} --pheno ../data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.${tissue}.final --mpheno $mpheno_number --thread-num 5 --out ../tem/${tissue}.${gene}
	h2=`grep 'V(G)/Vp' ../tem/${tissue}.${gene}.hsq`
	p=`grep 'Pval' ../tem/${tissue}.${gene}.hsq`
	if [ -n "${h2}" ] && [ -n "$p" ]; then
		echo -n $gene >> ../result/${tissue}.heritability
		echo -n ' ' >> ../result/${tissue}.heritability
		echo -n $h2 >> ../result/${tissue}.heritability
		echo -n ' ' >> ../result/${tissue}.heritability
		echo $p >> ../result/${tissue}.heritability
	fi
	rm ../tem/${tissue}.${gene}.hsq
	rm ../tem/${gene}.{bed,bim,fam,grm.bin,grm.id,grm.N.bin,log,nosex}
done
