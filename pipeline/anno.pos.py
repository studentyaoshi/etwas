import sys
a=open('../data/expression/gtex.gene.final.nomhc.anno','rt')
dic={}
for i in a:
	gene,chr,start,end,str=i.strip().split('\t')
	dic[gene]='\t'.join([gene,chr,start,end])
a.close()
a=open(sys.argv[1],'rt')
b=open(sys.argv[2],'wt')
b.write('WGT\tID\tCHR\tP0\tP1\n')
for i in a:
	gene=i.rstrip('\n')
	try:
		b.write(sys.argv[3]+'.rdata/'+gene+'.RData'+'\t'+dic.get(gene)+'\n')
	except:
		print gene
a.close()
b.close()
