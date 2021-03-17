tissue=$1 
tissuetype='brain'
genelist=../result/${tissue}.P005
prepheno='../data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
genegeno='../tem'

cat $genelist|while read gene
do 
	# Get best snp set
	geneperf=`grep $gene ../result/$tissue.etwas`
	bestmodel=`echo $geneperf|cut -d' ' -f 5`
	if [ $bestmodel = 'lasso' ];then
		snpset=`echo $geneperf|cut -d' ' -f 1`
	elif [ $bestmodel = 'enet' ];then
		snpset=`echo $geneperf|cut -d' ' -f 3`
	fi
	plink --bfile $genegeno/$tissue.$gene --assoc --pheno $genegeno/$tissue.$gene.expression --allow-no-sex --out ../tem/$tissue.$gene
	awk '$9<0.05{print$2"\t"$9}' ../tem/${tissue}.${gene}.qassoc > ../tem/${tissue}.${gene}.eqtl
	rm ../tem/$tissue.$gene.{qassoc,log,nosex}
	python getbestsnp.py $snpset $tissuetype $genegeno/$tissue.$gene.bim ../tem/$tissue.$gene.eqtl ../tem/$tissue.$gene.snp

	# Train the final models
	plink --bfile $genegeno/$tissue.$gene --extract ../tem/$tissue.$gene.snp --recodeA --recode-allele $genegeno/${tissue}.${gene}.minor --out ../tem/$tissue.$gene.train
	python Top_get.train.matrix.py $genegeno/${tissue}.${gene}.expression ../tem/$tissue.$gene.train.raw ../tem/$tissue.$gene.train
	
	rm ../tem/$tissue.$gene.train.{raw,nosex,log}
	if [ $bestmodel = 'lasso' ];then
		R CMD BATCH --slave "--args ../tem/$tissue.$gene.train 1 ../tem/$tissue.$gene.weight.txt" Get_snpweights.r
	elif [ $bestmodel = 'enet' ];then
		R CMD BATCH --slave "--args ../tem/$tissue.$gene.train 0.5 ../tem/$tissue.$gene.weight.txt" Get_snpweights.r
	fi
	
	python change_weight.py ../tem/$tissue.$gene.weight.txt ../tem/$tissue.$gene.weights
	cut -f 1 ../tem/$tissue.$gene.weights|while read snp
	do
		grep -w $snp $genegeno/$tissue.$gene.bim >> ../tem/$tissue.$gene.snps
	done

	if [ -s ../tem/$tissue.$gene.snps ];then
		R CMD BATCH --slave "--args ../tem/$tissue.$gene.weights ../tem/$tissue.$gene.snps ../result/$tissue.rdata/$gene.RData" get_rdata.r
	fi
	rm ../tem/$tissue.$gene.{eqtl,snp,snps,train,weight.txt}

	if [ -s ../result/$tissue.rdata/$gene.RData ];then
		echo $gene >> ../result/$tissue.finish.genes
	fi
done

## Get the rdata.path
python anno.pos.py ../result/$tissue.finish.genes ../result/$tissue.pos $tissue
