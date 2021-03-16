import sys
a=open(sys.argv[1],'rt')
dic={}
dic['IID']='expression'
a.readline()
for i in a:
	all=i.strip().split('\t')
	dic[all[1]]=all[2]
a.close()
a=open(sys.argv[2],'rt')
b=open(sys.argv[3],'wt')
for i in a:
	all=i.strip().split(' ',6)
	b.write(all[1]+','+all[6].replace(' ',',')+','+dic.get(all[1])+'\n')
a.close()
b.close()
