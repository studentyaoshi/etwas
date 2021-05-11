import sys
a=open(sys.argv[1],'rt')
zd={}
for i in a:
	all=i.strip().split('\t')
	s=all[1]
	zd[s]='\t'.join([all[0],all[1],all[2],all[3],all[4],all[5]])
a.close()
a=open(sys.argv[2],'rt')
for i in a:
	s=i.strip().split('\t')[1]
	print zd.get(s)
a.close()
