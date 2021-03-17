import sys
import random
a=open(sys.argv[1],'rt')
b=open(sys.argv[2],'wt')
a.readline()
a.readline()
for i in a:
	snp,w=i.strip().split(' ')
	if w!='0':
		b.write(snp.split('_')[0].split('"')[1]+'\t'+w+'\t'+str(random.uniform(-1,1))+'\n')
a.close()
b.close()
