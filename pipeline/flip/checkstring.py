import sys
a=open(sys.argv[1],'rt')
dic={}
dic['A']='T'
dic['T']='A'
dic['C']='G'
dic['G']='C'
x=range(1,23)
geno=[]
for w in x:
	filename='/home/ys/data/dbGaP/pipeline_impute_illu/genome/chr'+str(w)+'.fa'
	b=open(filename,'rt')
	geno.append(b.readline())
	b.close()
for i in a:
	try:
		all=i.strip().split(' ',5)
		chro=all[0].strip().split('r')[1]
		pos=all[1]
		pypos=int(int(pos)-1)
		chroo=int(chro)-1
		genol=geno[chroo]
		re=genol[pypos].upper()
		al1=all[4]
		al2=all[5]
		if al1!=al2 and al1!=dic.get(al2):
			if al1==re or al2==re:
				print '\t'.join([chro,all[3],'0',pos,all[4],all[5],'+'])
			elif al1!=re and al2!=re:
				print '\t'.join([chro,all[3],'0',pos,all[4],all[5],'-'])
	except:
#		print i,
		all=i.strip().split(' ',4)
		chro=all[0].strip().split('r')[1]
		pos=str(int(all[1])-1)
		pypos=int(int(pos)-1)
		chroo=int(chro)-1
		genol=geno[chroo]
		re=genol[pypos].upper()
		al1=all[3]
		al2=all[4]
		if al1!=al2 and al1!=dic.get(al2):
			if al1==re or al2==re:
				print '\t'.join([chro,all[2],'0',pos,all[3],all[4],'+'])
			elif al1!=re and al2!=re:
				print '\t'.join([chro,all[2],'0',pos,all[3],all[4],'-'])
a.close()
