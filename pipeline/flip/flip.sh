# Convert to bedtools format
python change_bim.py ref.bim > tem.bed

# Convert genome position to GRCh38
liftOver tem.bed hg19ToHg38.over.chain.gz output.bed unlifted.bed 
	## We provide hg19ToHg38 and hg18ToHg38 downloaded from http://hgdownload.soe.ucsc.edu/downloads.html. 

# Check position	
cat output.bed|while read LINE
do
	var=`echo "$LINE"|awk -F'\t' '{print$2}'`
	snp=`echo "$LINE"|awk -F'\t' '{print$4}'`
	let min=$var-100
	let max=$var+100
	chr=`echo "$LINE"|awk -F'\t' '{print$1}'|awk -F'r' '{print$2}'`
	pos=`tabix /home/ys/data/1000G_raw/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz $chr:$min-$max|grep ${snp}|awk -F'\t' '{print$2}'`
	if [ "$var" == "$pos" ] ;then
		echo $LINE >> output.change.bed
	else
		echo ${LINE/$var/$pos} >> output.change.bed
	fi
done
	## 

python checkstring.py output.change.bed > withchain.bim
awk -F '\t' {'print$2'} withchain.bim> snp.list
plink --bfile CIDR_AutopsyPD_TOP_subject_level --extract snp.list --make-bed --out snp_list
python /cu02_public/ys/IMPUTE/pipeline/pipeline_impute_illu/new_bim.py withchain.bim snp_list.bim > snp_list2.bim
mv snp_list2.bim snp_list.bim
plink --bfile snp_list --make-bed --out snp_list2
awk -F '\t' '$7=="-"{print$2}' withchain.bim > minus.snp
plink --bfile snp_list2 --flip minus.snp --make-bed --out flip_list
