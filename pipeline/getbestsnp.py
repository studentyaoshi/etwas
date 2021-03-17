import sys
snpset=sys.argv[1]
anno=snpset.strip().split('.',6)
hmm=anno[3]
dhs=anno[4]
tfbs=anno[5]
eqtl=anno[6]

hmmdic={}
a=open('../data/epigenetic/'+sys.argv[2]+'.'+hmm+'.snp','rt')
for i in a:
	hmmdic[i.rstrip('\n')]=1
a.close()

dhsdic={}
dhsyes={}
a=open('../data/epigenetic/'+sys.argv[2]+'.dnase.snp','rt')
for i in a:
	dhsyes[i.rstrip('\n')]=1
a.close()
if dhs=='dy':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		if dhsyes.get(rs)==1:
			dhsdic[rs]=1
	b.close()
elif dhs=='dn':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		if dhsyes.get(rs)!=1:
			dhsdic[rs]=1
	b.close()
elif dhs=='dydn':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		dhsdic[rs]=1
	b.close()
a.close()

tfbsdic={}
tfbsyes={}
a=open('../data/epigenetic/tfbs.snp','rt')
for i in a:
	tfbsyes[i.rstrip('\n')]=1
a.close()
if tfbs=='ty':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		if tfbsyes.get(rs)==1:
			tfbsdic[rs]=1
	b.close()
elif tfbs=='tn':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		if tfbsyes.get(rs)!=1:
			tfbsdic[rs]=1
	b.close()
elif tfbs=='tytn':
	b=open(sys.argv[3],'rt')
	for i in b:
		rs=i.strip().split('\t')[1]
		tfbsdic[rs]=1
	b.close()

eqtldic={}
a=open(sys.argv[4],'rt')
b=open(sys.argv[5],'wt')
for i in a:
	snp,p=i.strip().split('\t')
	if hmmdic.get(snp)==1 and dhsdic.get(snp)==1 and tfbsdic.get(snp)==1 and float(p)<float(eqtl):
		b.write(snp+'\n')
a.close()
b.close()
