tissue=$1
tissuetype='brain'
genelist=../result/${tissue}.P005
prepheno='../data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
genegeno='../tem'

cat $genelist|while read gene
do
	for i in {a..e}
	do
		# get eqtl
		plink --bfile $genegeno/$tissue.$gene --remove ${prepheno}.${tissue}.fam.a$i --assoc --pheno $genegeno/$tissue.$gene.expression --out ../tem/$tissue.$gene.$i --allow-no-sex
		awk '$9<0.05{print$2"\t"$9}' ../tem/${tissue}.${gene}.$i.qassoc > ../tem/${tissue}.${gene}.$i.eqtl
		rm ../tem/$tissue.$gene.$i.{qassoc,log,nosex}
		# get snp sets
		python getsnpset.py $tissuetype ../tem/$tissue.$gene.$i.eqtl ../tem/$tissue.$gene.$i
		rm ../tem/${tissue}.${gene}.$i.eqtl
		for eqtl in {'0.05','0.01','0.001','0.0001','0.00001','0.000001'}
		do 
			for hmm in {'enh','tss','tx','other','allhmm'}
			do 
				for tfbs in {'ty','tn','tytn'}
				do
					for dhs in {'dy','dn','dydn'}
					do 
						if [ -s ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.snp ]; then
							plink --bfile $genegeno/$tissue.$gene --remove ${prepheno}.${tissue}.fam.a$i --extract ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.snp --recodeA --recode-allele $genegeno/${tissue}.${gene}.minor --out ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.train
							python Top_get.train.matrix.py $genegeno/${tissue}.${gene}.expression ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.train.raw ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.train

							plink --bfile $genegeno/$tissue.$gene --keep ${prepheno}.${tissue}.fam.a$i --extract ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.snp --recodeA --recode-allele $genegeno/${tissue}.${gene}.minor --out ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.test
							python Top_get.train.matrix.py $genegeno/${tissue}.${gene}.expression ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.test.raw ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.test

							R CMD BATCH --slave "--args ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.train ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.$eqtl.test ${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}" Active.predixcan.r

							rm ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.{snp,test,test.log,test.nosex,test.raw,train,train.log,train.nosex,train.raw}
						else
							rm ../tem/${tissue}.${gene}.$i.${hmm}.${dhs}.${tfbs}.${eqtl}.snp
						fi
					done
				done
			done
		done
	done
	for eqtl in {'0.05','0.01','0.001','0.0001','0.00001','0.000001'}
	do
		for hmm in {'enh','tss','tx','other','allhmm'}
		do
			for tfbs in {'ty','tn','tytn'}
			do
				for dhs in {'dy','dn','dydn'}
				do
					if [ -e ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2 ]; then
						sed -i -e '/NA/d' ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2
						if [ -s ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2 ];then
							mean=`awk '{sum+=$1}END{mean=sum/5;printf "%8.8f\n",mean}' ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2`
							echo ${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}' '$mean >> ../tem/${tissue}.${gene}.lasso.r2
							rm ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2
						else
							rm ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.lasso.r2
						fi
					fi
					if [ -e ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2 ]; then
						sed -i -e '/NA/d' ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2
						if [ -s ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2 ];then
							mean=`awk '{sum+=$1}END{mean=sum/5;printf "%8.8f\n",mean}' ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2`
							echo ${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}' '$mean >> ../tem/${tissue}.${gene}.enet.r2
							rm ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2
						else
							rm ../tem/${tissue}.${gene}.${hmm}.${dhs}.${tfbs}.${eqtl}.enet.r2
						fi
					fi
				done
			done
		done
	done
	lassom=`sort -nrk2 ../tem/$tissue.$gene.lasso.r2|head -1|cut -f 1 -d' '`
	lassor=`sort -nrk2 ../tem/$tissue.$gene.lasso.r2|head -1|cut -f 2 -d' '`
	enetm=`sort -nrk2 ../tem/$tissue.$gene.enet.r2|head -1|cut -f 1 -d' '`
	enetr=`sort -nrk2 ../tem/$tissue.$gene.enet.r2|head -1|cut -f 2 -d' '`
	#rm ../tem/$tissue.$gene.{lasso,enet}.r2
	if [ `echo "$enetr > $lassor"|bc` -eq 1 ];then
		echo $lassom' '$lassor' '$enetm' '$enetr' enet '$enetr >> ../result/$tissue.etwas
	else
		echo $lassom' '$lassor' '$enetm' '$enetr' lasso '$lassor >> ../result/$tissue.etwas
	fi
done
