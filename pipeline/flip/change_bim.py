import sys
a=open(sys.argv[1],'rt')
for i in a:
	all=i.strip().split('\t')
	chr='chr'+all[0]
	sta=all[3]
	end=str(int(all[3])+1)
	if all[4]!='0' and all[5]!='0':
		print chr+'\t'+sta+'\t'+end+'\t'+all[1]+'\t'+all[4]+'\t'+all[5]
a.close()
