cat ../data/expression/gtex.gene.final.nomhc.anno|while read line
do gene=`echo $line|awk '{print$1}'`
	chr=`echo $line|awk '{print$2}'`
	position1=`echo $line|awk '{print$3}'`
	position2=`echo $line|awk '{print$4}'`
	let start=$position1-1000000
	let end=$position2+1000000
	if [ $start -lt 0 ]; then start=0; fi
	plink --bfile ../../genotype/GTEx.chr$chr --maf 0.05 --chr $chr --from-bp $start --to-bp $end --make-bed --out ../tem/$gene
	gcta64 --bfile ../tem/$gene --make-grm --out ../tem/${gene} --thread-num 5
	for tissue in {'Brain_Amygdala','Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia','Brain_Spinal_cord_cervical_c-1','Brain_Substantia_nigra','Whole_Blood'}
	do
		mpheno_number=`head -1 ../../expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.${tissue}.final|sed 's/\t/\n/g'|sed '1,2d' |grep -n -w ${gene}|awk -F ':' '{print$1}'`
		gcta64 --reml --grm ../tem/${gene} --pheno ../../expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.${tissue}.final --mpheno $mpheno_number --thread-num 5 --out ../tem/${tissue}.${gene}
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
	done
	rm ../tem/${gene}.{bed,bim,fam,grm.bin,grm.id,grm.N.bin,log,nosex}
done
