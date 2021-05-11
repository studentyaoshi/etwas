<div align=center>
<img src="./logo/bigc.png" width="400" height="200" slt="bigclogo" align="middle" />
</div>

# ETWAS
**E**pigenetic element-based **T**ranscriptome-**W**ide **A**ssociation **S**tudies
## Introduction  
> Transcriptome-wide association studies (TWAS) providing a powerful approach to identify novel disease risk genes and uncover possible causal genes at loci identified previously by genome-wide association study (GWAS). However, these methods did not consider the importance of epigenetic regulation in gene expression. Here, we developed a novel epigenetic element-based transcriptome-wide association study (ETWAS) that tested the effects of genetic variants on gene expression levels with the epigenetic features as prior and further mediated the association between predicted expression and complex diseases.  

## Contact
### Citation
> **Method**: **Yao Shi**, et al. [Epigenetic Element-Based Transcriptome-Wide Association Study Identifies Novel Genes for Bipolar Disorder](https://academic.oup.com/schizophreniabulletin/advance-article-abstract/doi/10.1093/schbul/sbab023/6190174?redirectedFrom=fulltext). *Schizophrenia Bulletin* 2021.  
> **Application**: **Yao Shi**, et al. A Transcriptome-wide Association Study Identifies Susceptibility Genes for Parkinson’s Disease. *Under Review*.  
> The data that support the findings of this study are available from the corresponding author upon reasonable request.
### Author
> **Yao Shi**, **Guo Yan**  
> National and Local Joint Engineering Research Center of Biodiagnosis and Biotherapy, The Second Affiliated Hospital, Xi'an Jiaotong University, Xi'an, Shaanxi, 710004, P. R. China  
> Key Laboratory of Biomedical Information Engineering of Ministry of Education, School of Life Science and Technology, Xi'an Jiaotong University, Xi'an, Shaanxi, 710049, P. R. China  
> :email:guoyan253@xjtu.edu.cn  
### Maintainer
> **Yao Shi**  
> You can contact :email:studentyaoshi@stu.xjtu.edu.cn when you have any questions, suggestions, comments, etc.  
> Please describe in details, and attach your command line and log messages if possible.  

## Requirements  
> - [**GCTA**](http://cnsgenomics.com/software/gcta/)
> - [**Bedtools**](http://quinlanlab.org/tutorials/bedtools/bedtools.html)
> - [**Python**](https://www.python.org/downloads/)
> - [**Plink**](http://zzz.bwh.harvard.edu/plink/epidetails.shtml)
> - [**FUSION**](http://gusevlab.org/projects/fusion/)
> - [**LDSC**](https://github.com/bulik/ldsc)
> - [**R**](https://www.r-project.org/)
>	- R packages: glmnet, doParallel.  

## Schematic of ETWAS
### Heritability calculation
Our current study is based on the premise that gene expression is heritable. Considering the heritability genes typically enriched for trait associations, we estimated gene expression heritability and only focused on the significantly heritable genes in further analyses.
- Need files
	- ~/etwas/data/genotype/GTEx.plink `genotype`
	
	The genotype in plink format can be found [**here**](http://zzz.bwh.harvard.edu/plink/data.shtml#bed).

	Note: All allele in this study were aligned on the forward strand and the SNP position were assigned to the GRCh38 human reference genome assembly. We provide the pipeline to flip the strand at [**here**](test).
	- ~/etwas/data/expression/$tissue.final `gene expression`
	```
	FID	IID	ENSG00000186092.4	ENSG00000278566.1	ENSG00000273547.1	ENSG00000187634.11	ENSG00000188976.10	ENSG00000187961.13	...
	GTEX-13O3Q	GTEX-13O3Q	0.2272	0.1388	0.02777	0.49	30.3	4.328	...
	GTEX-1C6VR	GTEX-1C6VR	0.1296	0.1267	0.03167	0.9371	13.5	2.371	...
	GTEX-P44G	GTEX-P44G	0.3325	0.266	0.1478	0.7702	19.54	3.004	...
	GTEX-X4XX	GTEX-X4XX	0.07441	0.1455	0.1455	1.323	26.77	11.23	...
	GTEX-1EX96	GTEX-1EX96	0.07624	0.04969	0	0.3575	8.479	1.611	...
	```
	- ~/etwas/data/expression/gtex.gene.final.nomhc.anno `gene annotation`
	```
	ENSG00000186092.4	1	69091	70008	+
	ENSG00000278566.1	1	450740	451678	-
	ENSG00000273547.1	1	685716	686654	-
	ENSG00000187634.11	1	923928	944581	+
	ENSG00000188976.10	1	944204	959309	-
	```
- Run
```
cd ~/etwas/pipeline
sh heritability.sh $tissue

awk '$6<0.05{print$1}' ~/etwas/result/$tissue.heritability > ~/etwas/result/$tissue.P005
```
The $tissue indicates the tissue name, such as `Brain_Amygdala`.
- Generate
	- ~/etwas/result/$tissue.heritability
	```
	ENSG00000177989.13 V(G)/Vp 0.166964 0.113988 Pval 0.03208
	ENSG00000130487.5 V(G)/Vp 0.063453 0.091522 Pval 0.218
	ENSG00000217442.3 V(G)/Vp 0.177535 0.114586 Pval 0.02242
	ENSG00000205560.12 V(G)/Vp 0.049175 0.082262 Pval 0.233
	ENSG00000100288.19 V(G)/Vp 0.000001 0.063587 Pval 0.5
	```
	- ~/etwas/result/$tissue.P005 `gene list with heritability p-value less than 0.05`  

We estimated the cis heritability (1 Mb window around each gene) for each gene using restricted maximum likelihood analysis, a variance-component model with a genetic relationship matrix (GRM) estimated from genotype data in GCTA software. Genes with heritability *P*-value less than 0.05 were regarded as significantly heritable genes.

### Model generation
For each gene *x*, we first included SNPs within the 1 MB region around gene. We trained and evaluated the models for gene expression prediction in each round *y* of tenfold cross-validation using the following steps:  
1. We performed eQTL analyses with SNPs located within 1 Mb of the transcription start/end sites of the gene using the training data.  
2. We then annotated the SNPs with epigenetic annotations. The epigenomic data included chromatin segmentation states, transcription factor binding sites (TFBS), and DNase I hypersensitive sites (DHS). For the chromatin segmentation states, we utilized the 15-state model and grouped them into four categories: *promoter*, *enhancer*, *transcription*, and *others*. Individually, *TssA* and *TssAFlnk* were considered to be a *promoter*; *TxFlnk*, *Tx*, and *TxWk* were considered to be a *transcription*; *EnhG* and *Enh* were considered to be an *enhancer*, and the rest were classified as the *others*.  
	For each SNP, an epigenomic feature was labeled if the SNP overlapped with the feature.  
3. We obtained multiple SNP sets according to the eQTL *P*-value threshold (5×10<sup>-2</sup>, 1×10<sup>-2</sup>, 1×10<sup>-3</sup>, 1×10<sup>-4</sup>, 1×10<sup>-5</sup>, 1×10<sup>-6</sup>), chromatin segmentation state (*promoter*, *transcription*, *enhancer*, and *others*), TFBS annotation, and DHS annotation.  
4. For each SNP set *z*, we built an expression prediction model in the training dataset by using the lasso and the elastic net (α = 0.5) methods as implemented in the *glmnet* R package.
- Need files
	- ~/etwas/result/$tissue.P005 `gene list with heritability p-value less than 0.05`
	- ~/etwas/data/genotype/GTEx.plink `genotype`
	- ~/etwas/data/expression/$tissue.final `gene expression`
	- ~/etwas/data/expression/gtex.gene.final.nomhc.anno `gene annotation`
	- ~/etwas/data/epigenetic/\* `epigenetic annotation`
	```
	chr1	10131	10369
	chr1	10427	10574
	chr1	57344	57408
	chr1	235654	235786
	chr1	235887	235971
	```
- Run
```
# Get genotype for each gene
cd ~/etwas/pipeline
sh genopheno.sh $tissue

# Epigenetic-based SNP grouping and model training
sh train.sh $tissue
```
- Generate
	- ~/etwas/tem/$tissue.$gene.plink `genotype for each gene`
	- ~/etwas/tem/$tissue.$gene.expression `expression for each gene`
	```
	> Brain_Amygdala.ENSG00000024526.16.expression
	FID	IID	ENSG00000024526.16
	GTEX-13O3Q	GTEX-13O3Q	0.03047
	GTEX-1C6VR	GTEX-1C6VR	0.02316
	GTEX-P44G	GTEX-P44G	0.08647
	GTEX-X4XX	GTEX-X4XX	0.02661
	GTEX-1EX96	GTEX-1EX96	0
	```
	- ~/etwas/tem/$tissue.$gene.minor `minor alleles for SNPs in gene`
	
	Note: We assigned the minor allele (in all reference samples) to be the reference allele (A1) to avoid the allele difference among cross-validation.
	```
	> Brain_Amygdala.ENSG00000024526.16.minor
	SNP	A1
	rs1977785	G
	rs12745192	G
	rs6695089	A
	rs7554654	G
	rs11209084	G
	```
	- ~/etwas/result/$tissue.etwas
	```
	Brain_Amygdala.ENSG00000169885.9.enh.dn.tn.0.05 0.09449683 Brain_Amygdala.ENSG00000169885.9.enh.dn.tn.0.05 0.09388865 lasso 0.09449683
	Brain_Amygdala.ENSG00000116151.13.tx.dn.tn.0.05 0.12374554 Brain_Amygdala.ENSG00000116151.13.tx.dn.tn.0.05 0.12354224 lasso 0.12374554
	Brain_Amygdala.ENSG00000157916.19.tx.dn.tn.0.05 0.11186665 Brain_Amygdala.ENSG00000157916.19.tx.dn.tn.0.05 0.11255954 enet 0.11255954
	Brain_Amygdala.ENSG00000157881.13.enh.dn.tn.0.05 0.12066406 Brain_Amygdala.ENSG00000157881.13.enh.dn.tn.0.05 0.12117371 enet 0.12117371
	Brain_Amygdala.ENSG00000157873.17.tx.dn.tn.0.05 0.17463774 Brain_Amygdala.ENSG00000157873.17.tx.dn.tn.0.05 0.23938381 enet 0.23938381
	```
### Get final model
For each model, we evaluated its prediction performance by cross-validation R<sup>2</sup> between the predicted gene expression and the observed gene expression of the testing data and averaged all the cross-validation data. For each gene *x*, the model with the highest mean R<sup>2</sup> in the testing data was selected as the best model. Based on the parameters of the best model, we performed the eQTL analyses again using all the samples in the reference data and constructed each gene’s final prediction model.
- Need files:
	- ~/etwas/result/$tissue.P005 `gene list with heritability p-value less than 0.05`
	- ~/etwas/result/$tissue.etwas `etwas results`
	- ~/etwas/tem/$tissue.$gene.plink `genotype for each gene`
	- ~/etwas/tem/$tissue.$gene.expression `expression for each gene`
- Run
```
cd ~/etwas/pipeline
sh getfinal.sh $tissue
```
- Generate
	- ~/etwas/result/$tissue.pos
	```
	WGT	ID	CHR	P0	P1	N
	Brain_Amygdala.rdata/ENSG00000169885.9.RData	ENSG00000169885.9   1	1914827	1917296	120
	Brain_Amygdala.rdata/ENSG00000116151.13.RData	ENSG00000116151.13  1	2321253	2391707	120
	Brain_Amygdala.rdata/ENSG00000157916.19.RData	ENSG00000157916.19  1	2391775	2405444	120
	Brain_Amygdala.rdata/ENSG00000157881.13.RData	ENSG00000157881.13  1	2508533	2526628	120
	Brain_Amygdala.rdata/ENSG00000215912.12.RData	ENSG00000215912.12  1	2635976	2801717	120
	```
	- ~/etwas/result/$tissue.finish.genes
	```
	ENSG00000001167.14
	ENSG00000004975.11
	ENSG00000005339.14
	ENSG00000005469.11
	ENSG00000005471.15
	```
	- ~/etwas/result/$tissue.rdata/ `this folder includes all the .RData files in the $tissue`
	```
	ENSG00000001167.14.RData
	ENSG00000004975.11.RData
	ENSG00000005339.14.RData
	ENSG00000005469.11.RData
	ENSG00000005471.15.RData
	```
	each $gene.RData includes "wgt.matrix" and "snps"
	```
	> wgt.matrix
			etwas       top1
	rs2916256   2.46794091 -0.3117201
	rs9394715  -1.15137795 -0.1625038
	rs12664653 -0.34564314 -0.2289465
	rs10947907 -0.26958668  0.1971410
	rs4521573  -0.11433455 -0.7497704
	rs4236071  -0.12000062 -0.7489674
	rs9394733  -0.06642380  0.7040359
	rs7767835  -0.02499029  0.4008386
	rs9367079   1.38613939  0.1223406
	> snps
	  V1         V2 V3       V4 V5 V6
	1  6  rs2916256  0 40383313  C  T
	2  6  rs9394715  0 40659686  C  A
	3  6 rs12664653  0 40673177  G  A
	4  6 rs10947907  0 40674904  C  T
	5  6  rs4521573  0 40773205  C  T
	6  6  rs4236071  0 40782300  T  C
	7  6  rs9394733  0 40785517  C  T
	8  6  rs7767835  0 40793076  G  T
	9  6  rs9367079  0 40808854  G  A
	```
We provided all the needed files and results in the **Brain_Amygdala**, the results in other Brain tissues will be uploaded soon.
### TWAS
After getting the best model for gene *x*, we could predict expression directly for genotyped samples using the effect sizes from the reference panels and measure the association between predicted expression and a trait. On the other hand, the [ImpG-Summary](https://academic.oup.com/bioinformatics/article/30/20/2906/2422225) algorithm has been used to extend to train on the genetic component of expression based on GWAS summary data. Thus, we could indirectly estimate the association between predicted expression and a trait as the weighted linear combination of SNP-trait standardized effect sizes while accounting for linkage disequilibrium (LD) among SNPs. [FUSION](http://gusevlab.org/projects/fusion/) was used to conduct the transcriptome-wide association testing. The 1000 Genomes v3 LD panel was used for the ETWAS.
- Need files
	- ~/etwas/result/$tissue.pos `etwas results`
	- ~/etwas/result/$tissue.rdata/ `etwas results`
	- ~/etwas/data/LDREF/ `ld reference data`
	- ~/etwas/data/gwas/$trait.sumstats.gz `summary statistic of trait`
	
	The summary statistic file format can be found [**here**](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format).
- Run
```
sh predict.sh $tissue $trait 
```
The $trait indicates the name of the summary data in ~/etwas/data/gwas/, such as `PGC_BIP2018`.
- Generate
	- ~/etwas/result/twas/$trait.$tissue.$chr
	```
	> PGC_BIP2018.Brain_Amygdala.1
	PANEL	FILE	ID	CHR	P0	P1	HSQ	BEST.GWAS.ID	BEST.GWAS.Z	EQTL.ID	EQTL.R2	EQTL.Z	EQTL.GWAS.Z	NSNP	NWGT	MODEL	MODELCV.R2    MODELCV.PV	TWAS.Z	TWAS.P	PERM.PV	PERM.N	PERM.ANL_PV
	Brain_Amygdala	~/etwas/result/Brain_Amygdala.rdata/ENSG00000169885.9.RData ENSG00000169885.9	1	1.91e+06	1.92e+06	0	rs34238514	-1.847	rs12041925	0	-0.8795	 0.588   7	 7		0	0	 1.31872	0.18726	0.0000	    0	 0.00e+00
	Brain_Amygdala	~/etwas/result/Brain_Amygdala.rdata/ENSG00000116151.13.RData	ENSG00000116151.13	1	2.32e+06	2.39e+06	0	rs2500271	-2.840	rs4531246	0	 0.9934	-0.306	13	13		0	0	-0.42751	0.66901	0.0000	    0	 0.00e+00
	Brain_Amygdala	~/etwas/result/Brain_Amygdala.rdata/ENSG00000157916.19.RData	ENSG00000157916.19	1	2.39e+06	2.41e+06	0	rs1108600	-2.051	rs12410859	0	 0.9923	-1.707	16	16		0	0	-0.70353	0.48172	0.0000	    0	 0.00e+00
	Brain_Amygdala	~/etwas/result/Brain_Amygdala.rdata/ENSG00000157881.13.RData	ENSG00000157881.13	1	2.51e+06	2.53e+06	0	rs1108600	-2.051	rs6603813	0	-0.9847	-1.729	 9	 9		0	0	-1.66058	0.09680	0.0000	    0	 0.00e+00
	Brain_Amygdala	~/etwas/result/Brain_Amygdala.rdata/ENSG00000215912.12.RData	ENSG00000215912.12	1	2.64e+06	2.80e+06	0	rs2503701	 1.494	rs2503701	0	 0.2924	 1.494	 2	 2		0	0	-0.49232	0.62250	0.0000	    0	 0.00e+00
	```
