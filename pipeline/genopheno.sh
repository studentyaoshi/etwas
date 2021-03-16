tissue=$1 
genelist=../result/${tissue}.P005
prepheno='../data/expression/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct'
pregeno='../data/genotype'
cat $genelist|while read gene
do
	# Get the basic information of gene
	anno=`grep -w $gene ../data/expression/gtex.gene.final.nomhc.anno`
	chr=`echo $anno|cut -d' ' -f 2`
	tss=`echo $anno|cut -d' ' -f 3`
	tes=`echo $anno|cut -d' ' -f 4`
	let start=$tss-1000000
	let end=$tes+1000000
	if [ $start -lt 0 ];then
		start=0
	fi

	#get expression
	number=`head -1 ${prepheno}.${tissue}.final|awk -v var=${gene} '{for(i=1;i<=NF;i++)if($i==var)print i}'`
	cut -f 1,2,$number ${prepheno}.${tissue}.final > ../tem/${tissue}.${gene}.expression
	#get genotype
	plink --bfile ${pregeno}/GTEx.chr$chr --chr $chr --from-bp $start --to-bp $end --allow-no-sex --keep ${prepheno}.${tissue}.fam --make-bed --out ../tem/${tissue}.${gene}
	plink --bfile ../tem/${tissue}.${gene} --freq --out ../tem/${tissue}.${gene}.freq
	awk '{print $2"\t"$3}' ../tem/${tissue}.${gene}.freq.frq > ../tem/${tissue}.${gene}.minor
	rm ../tem/$tissue.$gene.{freq.frq,freq.log,freq.nosex,log,nosex}
done
