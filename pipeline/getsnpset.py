import sys
eqtlthreshold=['0.05','0.01','0.001','0.0001','0.00001','0.000001']
hmmlist=['enh','tss','tx','other','allhmm']
tfbslist=['ty','tn',['ty','tn']]
dhslist=['dy','dn',['dy','dn']]

# hmm
enhdic={}
a=open('../data/epigenetic/'+sys.argv[1]+'.enh.snp','rt')
for i in a:
	enhdic[i.rstrip('\n')]=1
a.close()
tssdic={}
a=open('../data/epigenetic/'+sys.argv[1]+'.tss.snp','rt')
for i in a:
	tssdic[i.rstrip('\n')]=1
a.close()
txdic={}
a=open('../data/epigenetic/'+sys.argv[1]+'.tx.snp','rt')
for i in a:
	txdic[i.rstrip('\n')]=1
a.close()
otherdic={}
a=open('../data/epigenetic/'+sys.argv[1]+'.other.snp','rt')
for i in a:
	otherdic[i.rstrip('\n')]=1
a.close()
hmmdic={}
a=open(sys.argv[2],'rt')
for i in a:
	snp=i.strip().split('\t')[0]
	enh=enhdic.get(snp)
	tss=tssdic.get(snp)
	tx=txdic.get(snp)
	other=otherdic.get(snp)
	x=['allhmm']
	if enh==1:
		x.append('enh')
	if tss==1:
		x.append('tss')
	if tx==1:
		x.append('tx')
	if other==1:
		x.append('other')
	hmmdic[snp]=x

#dhs
dns1={}
a=open('../data/epigenetic/'+sys.argv[1]+'.dnase.snp','rt')
for i in a:
	dns1[i.rstrip('\n')]='1'
a.close()
dhsdic={}
a=open(sys.argv[2],'rt')
for i in a:
	snp=i.strip().split('\t')[0]
	if dns1.get(snp)=='1':
		dhsdic[snp]='dy'
	else:
		dhsdic[snp]='dn'
a.close()

#tfbs
tfbs1={}
a=open('../data/epigenetic/tfbs.snp','rt')
for i in a:
	tfbs1[i.rstrip('\n')]='1'
a.close()
tfbsdic={}
a=open(sys.argv[2],'rt')
for i in a:
	snp=i.strip().split('\t')[0]
	if tfbs1.get(snp)=='1':
		tfbsdic[snp]='ty'
	else:
		tfbsdic[snp]='tn'
a.close()
#genesets
for eqtl in eqtlthreshold:
	for hmm in hmmlist:
		for dhs in dhslist:
			dhsname=''.join(str(e) for e in dhs)
			for tfbs in tfbslist:
				tfbsname=''.join(str(y) for y in tfbs)
				a=open(sys.argv[2],'rt')
				b=open(sys.argv[3]+'.'+hmm+'.'+dhsname+'.'+tfbsname+'.'+eqtl+'.snp','wt')
				for x in a:
					snp,p=x.strip().split('\t')
					if hmm in hmmdic.get(snp) and float(p)<float(eqtl) and tfbsdic.get(snp) in tfbs and dhsdic.get(snp) in dhs:
						b.write(snp+'\n')
				a.close()
				b.close()
