import sys
import subprocess
from sequence_lib import read_fasta

#ref = ['ACGT','TAGC','GGTC','ACCG','TGGA','ATTG','AGTC','CTGT']
#est = ['AAAA','CCCC','GGGG','TTTT','TGCA','ACGT','ACTG','CAGT']

#rname,ref = read_fasta("out")
#ename,est = read_fasta("out1")

#The Fastsp function
def FastSP(ref,est):
	name_list= []
	n = len(ref)
	m = len(est)
	if n != m:
		print("Input lists have different sizes")
	f = open('ref.aln', 'w')
	for i in range(n):
		name_list.append('>'+'L'+ str(i)+'0273514')
#	print name_list
	for i in range(n):
		f.write(name_list[i])
		f.write('\n')
		f.write(ref[i])
		f.write('\n')
	f.close()
	f = open('est.aln', 'w')
	for i in range(n):
                f.write(name_list[i])
                f.write('\n')
                f.write(est[i])
                f.write('\n')
	f.close()
	s = subprocess.check_output(["java","-jar","FastSP/FastSP.jar","-r","ref.aln","-e","est.aln"],stderr=subprocess.STDOUT).split()
	#print(s)
	return int(s[24]), int(s[32]), int(s[40])

#shared,ref,est = FastSP(ref,est)
#print(shared,ref,est)
#print(Modeler)
