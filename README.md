<div align=center>
<img src="./logo/bigc.png" width="400" height="200" slt="bigclogo" align="middle" />
</div>

# ETWAS
Epigenetic element-based Transcriptome-Wide Association Studies
## Introduction  
> Transcriptome-wide association studies (TWAS) providing a powerful approach to identify novel disease risk genes and uncover possible causal genes at loci identified previously by genome-wide association study (GWAS). However, these methods did not consider the importance of epigenetic regulation in gene expression. Here, we developed a novel epigenetic element-based transcriptome-wide association study (ETWAS) that tested the effects of genetic variants on gene expression levels with the epigenetic features as prior and further mediated the association between predicted expression and complex diseases.  

## Contact
### Citation
> Yao Shi, et al. [Epigenetic Element-Based Transcriptome-Wide Association Study Identifies Novel Genes for Bipolar Disorder](https://www.medrxiv.org/content/10.1101/2020.07.23.20161174v3). medRxiv 2020.  
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
> - [**R**](https://www.r-project.org/)
>	- R packages: glmnet, doParallel.  

## Schematic of ETWAS
### Heritability calculation
Our current study is based on the premise that gene expression is heritable. Considering the heritability genes typically enriched for trait associations, we estimated gene expression heritability and only focused on the significantly heritable genes in further analyses.
- Need files:
	- ~/etwas/data/genotype/GTEx.plink `genotype`
	- ~/etwas/data/expression/$tissue.final `gene expression`
	- ~/etwas/data/expression/gtex.gene.final.nomhc.anno `gene annotation`
- Run
```
cd ~/etwas/pipeline
sh heritability.sh $tissue
```
The $tissue indicates the tissue name, such as Brain_Amygdala.
- Generate
	- ~/etwas/result/$tissue.heritability
```
ENSG00000177989.13 V(G)/Vp 0.166964 0.113988 Pval 0.03208
ENSG00000130487.5 V(G)/Vp 0.063453 0.091522 Pval 0.218
ENSG00000217442.3 V(G)/Vp 0.177535 0.114586 Pval 0.02242
ENSG00000205560.12 V(G)/Vp 0.049175 0.082262 Pval 0.233
ENSG00000100288.19 V(G)/Vp 0.000001 0.063587 Pval 0.5
```
We estimated the cis heritability (1 Mb window around each gene) for each gene using restricted maximum likelihood analysis, a variance-component model with a genetic relationship matrix (GRM) estimated from genotype data in GCTA software. Genes with heritability P-value less than 0.05 were regarded as significantly heritable genes.

### Train gene expression model
For each gene *x*, we first included SNPs within the 1 MB region around gene. We trained and evaluated the models for gene expression prediction in each round *y* of tenfold cross-validation using the following steps:  
1. We performed eQTL analyses with SNPs located within 1 Mb of the transcription start/end sites of the gene using the training data.  
2. We then annotated the SNPs with epigenetic annotations. The epigenomic data included chromatin segmentation states, transcription factor binding sites (TFBS), and DNase I hypersensitive sites (DHS). For the chromatin segmentation states, we utilized the 15-state model and grouped them into four categories: *promoter*, *enhancer*, *transcription*, and *others*. Individually, *TssA* and *TssAFlnk* were considered to be a *promoter*; *TxFlnk*, *Tx*, and *TxWk* were considered to be a *transcription*; *EnhG* and *Enh* were considered to be an *enhancer*, and the rest were classified as the *others*.  
	For each SNP, an epigenomic feature was labeled if the SNP overlapped with the feature.  
3. We obtained multiple SNP sets according to the eQTL *P*-value threshold (5×10<sup>-2</sup>, 1×10<sup>-2</sup>, 1×10<sup>-3</sup>, 1×10<sup>-4</sup>, 1×10<sup>-5</sup>, 1×10<sup>-6</sup>), chromatin segmentation state (*promoter*, *transcription*, *enhancer*, and *others*), TFBS annotation, and DHS annotation.  
4. For each SNP set *z*, we built an expression prediction model in the training dataset by using the lasso and the elastic net (α = 0.5) methods as implemented in the *glmnet* R package.
- Need files:
	- ~/etwas/result/${tissue}.P005 `gene list with heritability p-value less than 0.05`
	- ~/etwas/data/genotype/GTEx.plink `genotype`
	- ~/etwas/data/expression/$tissue.final `gene expression`
	- ~/etwas/data/expression/gtex.gene.final.nomhc.anno `gene annotation`
- Run
```
# Get genotype for each gene
cd ~/etwas/pipeline
sh genopheno.sh $tissue

# Epigenetic-based SNP grouping and model training
sh train.sh $tissue
```
- Generate

